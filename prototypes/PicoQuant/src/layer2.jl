import Random.shuffle
import Base: push!

export random_contraction_plan
export contraction_plan_to_json, contraction_plan_from_json
export contract_pair!, contract_network!
export AbstractExecuter, DSLWriter, InteractiveExecuter, push!
export compress_tensor_chain!

using HDF5

"""
    function random_contraction_plan(network::TensorNetworkCircuit)

Function to create a random contraction plan
"""
function random_contraction_plan(network::TensorNetworkCircuit)
    closed_edges = [x for (x, y) in pairs(edges(network)) if y.src != nothing
                                                          && y.dst != nothing]
    shuffle(closed_edges)
end

# function cost_flops(network, plan)
# end
#
# function cost_max_memory(network, plan)
# end
#
# function find_contraction_plan(network, cost_function)
# end


"""
    function contraction_plan_to_json(plan::Array)

Function to serialise the contraction plan to json format
"""
function contraction_plan_to_json(plan::Array{Symbol})
    JSON.json(plan)
end

"""
    function contraction_plan_from_json(str::String)

Function to deserialize the contraction plan from a json string
"""
function contraction_plan_from_json(str::String)
    [Symbol(x) for x in JSON.parse(str)]
end

# *************************************************************************** #
#                             Executer types
# *************************************************************************** #

"Abstract type defining executers which determine the behavior of contract
functions"
abstract type AbstractExecuter end

"The dsl writer will get the contract functions to write dsl commands."
mutable struct DSLWriter <: AbstractExecuter
    dsl_filename::String
    tensor_data_filename::String
    output_data_filename::String

    function DSLWriter(dsl::String="contract_network.tl",
                       tensor_data="tensor_data.h5"::String,
                       output::String="")
        # create empty dsl file
        close(open(dsl, "w"))
        new(dsl, tensor_data, output)
    end
end

function push!(dsl::DSLWriter, instruction::String)
    open(dsl.dsl_filename, "a") do io
        write(io, instruction * "\n")
    end
end

"Interactive executer which gets the contract functions to act on tensor data
in a TensorNetworkCircuit interactively."
mutable struct InteractiveExecuter <: AbstractExecuter
    use_gpu::Bool
    memory_size_mb::Int64

    InteractiveExecuter(use_gpu::Bool=false,
                        memory_size_mb::Int64=0) = new(use_gpu, memory_size_mb)
end

# *************************************************************************** #
#                  Functions to contract a tensor network
# *************************************************************************** #

"""
    function contract_indices(A::Node, B::Node)

This function divides all of the indices of two nodes A and B into two arrays,
one array for all the indices that are shared between A and B and the second
array will have all the indices that are not shared.
"""
function contract_indices(A::Node, B::Node)
    # Create an array of index labels for the new node
    common_indices = intersect(A.indices, B.indices)
    remaining_indices_A = setdiff(A.indices, common_indices)
    remaining_indices_B = setdiff(B.indices, common_indices)
    uncommon_indices = union(remaining_indices_A, remaining_indices_B)

    return common_indices, uncommon_indices
end

"""
    function create_ncon_indices(A::Node, B::Node,
                                 common_indices::Array{Symbol, 1},
                                 uncommon_indices::Array{Symbol, 1})

This function returns ncon indices for the nodes A and B. These can then to be
passed to ncon in order to contract the tensor data associated with A and B.
"""
function create_ncon_indices(A::Node, B::Node,
                             common_indices::Array{Symbol, 1},
                             uncommon_indices::Array{Symbol, 1})

    common_index_map = Dict([(x[2], x[1]) for x in enumerate(common_indices)])
    uncommon_indices_map = Dict([(x[2], -x[1]) for x
                                    in enumerate(uncommon_indices)])

    index_map = merge(common_index_map, uncommon_indices_map)
    A_ncon_indices = [index_map[x] for x in A.indices]
    B_ncon_indices = [index_map[x] for x in B.indices]

    return A_ncon_indices, B_ncon_indices
end

"""
    function contract_network!(network::TensorNetworkCircuit,
                               plan::Array{Symbol, 1},
                               executer::DSLWriter,
                               output_shape::Union{String, Array{<:Integer, 1}})

Function to write dsl commands describing how to contract the given network
according to the given contraction plan.
"""
function contract_network!(network::TensorNetworkCircuit,
                           plan::Array{Symbol, 1},
                           executer::DSLWriter,
                           output_shape::Union{String, Array{<:Integer, 1}}="")

    tensor_data_filename = executer.tensor_data_filename
    output_data_filename = executer.output_data_filename
    if output_data_filename == ""
        output_data_filename = tensor_data_filename
    end

    # Write the tensor data contained in the tensor network to a file.
    # Append load commands to the dsl script for each tensor.
    h5open(tensor_data_filename, "w") do file
        for (node_label, node) in pairs(network.nodes)
            node_name = string(node_label)
            write(file, node_name, node.data)
            push!(executer, "tensor $node_name")
        end
    end

    # Loop through the plan and contract each edge in sequence
    for edge in plan
        # sometimes edges get contracted ahead of time if connecting two
        # tensors being contracted
        if edge in keys(network.edges)
            contract_pair!(network, edge, executer)
        end
    end

    # Contract disjoint pieces of the network, if any
    while length(network.nodes) > 1
        n1, state = iterate(network.nodes)
        n2, _ = iterate(network.nodes, state)
        contract_pair!(network, n1.first, n2.first, executer)
    end

    # Reshape the final tensor if a shape is specified by the user.
    if output_shape == "vector"
        vector_length = 2^length(network.output_qubits)
        push!(executer, "reshape node_$(network.counters["node"]) $(vector_length)")

    elseif output_shape != ""
        push!(executer, "reshape node_$(network.counters["node"]) " * join(output_shape, ","))
    end

    # Add a command to save the contracted the network under the
    # group name 'result'.
    push!(executer, "save node_$(network.counters["node"]) $output_data_filename result")
end

"""
    function contract_network!(network::TensorNetworkCircuit,
                               plan::Array{Symbol, 1},
                               executer::InteractiveExecuter,
                               output_shape::Union{String, Array{<:Integer, 1}})

Function to interactively contract the given network according to the given
contraction plan.
"""
function contract_network!(network::TensorNetworkCircuit,
                           plan::Array{Symbol, 1},
                           executer::InteractiveExecuter,
                           output_shape::Union{String, Array{<:Integer, 1}}="")

    # Loop through the plan and contract each edge in sequence.
    for edge in plan
        # Sometimes edges get contracted ahead of time if connecting two
        # tensors being contracted.
        if edge in keys(network.edges)
            contract_pair!(network, edge, executer)
        end
    end

    # Contract disjoint pieces of the network, if any.
    while length(network.nodes) > 1
        n1, state = iterate(network.nodes)
        n2, _ = iterate(network.nodes, state)
        contract_pair!(network, n1.first, n2.first, executer)
    end

    # Get the node corresponding to the fully contracted tensor network.
    node, _ = iterate(values(network.nodes))

    # Reshape the final tensor if a shape is specified by the user.
    if output_shape == "vector"
        vector_length = 2^length(network.output_qubits)
        node.data[:] = reshape_tensor(node.data, vector)

    elseif output_shape != ""
        node.data[:] = reshape_tensor(node.data, output_shape)
    end
end

"""
    function contract_pair!(network::TensorNetworkCircuit,
                            edge::Symbol,
                            executer::AbstractExecuter)

Contract a pair of nodes connected by the given edge.
"""
function contract_pair!(network::TensorNetworkCircuit,
                        edge::Symbol,
                        executer::AbstractExecuter)
    if edge in keys(network.edges)
        e = network.edges[edge]
        return contract_pair!(network, e.src, e.dst, executer)
    end
end

"""
    function contract_pair!(network::TensorNetworkCircuit,
                            A_label::Symbol,
                            B_label::Symbol,
                            executer::DSLWriter)

Function to create dsl commands which contract the nodes A and B of the network.
"""
function contract_pair!(network::TensorNetworkCircuit,
                        A_label::Symbol,
                        B_label::Symbol,
                        executer::DSLWriter)

    # Should only contract different tensors, so skip contraction if A = B.
    if A_label == B_label
        return nothing
    end

    # Get the tensors being contracted and denote them A and B
    A = network.nodes[A_label]
    B = network.nodes[B_label]

    # Sort the indices of A and B into shared and unshared.
    common_indices, remaining_indices = contract_indices(A, B)

    # Create array of ncon indices for A and B.
    A_ncon_indices, B_ncon_indices = create_ncon_indices(A, B,
                                                         common_indices,
                                                         remaining_indices)

    # Create and add the new node to the tensor network.
    C_label = new_label!(network, "node")
    C = Node(remaining_indices, zeros(0))
    network.nodes[C_label] = C

    # remove contracted edges
    for index in common_indices
        delete!(network.edges, index)
    end

    # replumb existing edges to point to new node
    # TODO: might become time consuming for large networks
    for index in remaining_indices
        e = network.edges[index]
        if e.src in [A_label, B_label]
            e.src = C_label
        elseif e.dst in [A_label, B_label]
            e.dst = C_label
        end
    end

    # delete the nodes that were contracted
    delete!(network.nodes, A_label)
    delete!(network.nodes, B_label)

    # Create and return the dsl commands for performing this contraction and
    # deleting the old tensors.
    ncon_str = "ncon " * string(C_label) * " "
    ncon_str *= string(A_label) * " " * join(A_ncon_indices, ",") * " "
    ncon_str *= string(B_label) * " " * join(B_ncon_indices, ",")
    push!(executer, ncon_str)

    push!(executer, "del " * string(A_label))
    push!(executer, "del " * string(B_label))
    C_label
end

"""
    function contract_pair!(network::TensorNetworkCircuit,
                            A_label::Symbol,
                            B_label::Symbol,
                            executer::InteractiveExecuter)

Function to interactively contract a pair of nodes A and B.
"""
function contract_pair!(network::TensorNetworkCircuit,
                        A_label::Symbol,
                        B_label::Symbol,
                        executer::InteractiveExecuter)

    # Should only contract different tensors, so skip contraction if A = B.
    if A_label == B_label
        return nothing
    end

    # Get the tensors being contracted and denote them A and B
    A = network.nodes[A_label]
    B = network.nodes[B_label]

    # Sort the indices of A and B into shared and unshared.
    common_indices, remaining_indices = contract_indices(A, B)

    # Create array of ncon indices for A and B.
    A_ncon_indices, B_ncon_indices = create_ncon_indices(A, B,
                                                         common_indices,
                                                         remaining_indices)

    # Get the tensor data and create a new label for the tensor created
    # by contracting A and B.
    C_label = new_label!(network, "node")
    C_data = contract_tensors((A.data, B.data),
                              (A_ncon_indices, B_ncon_indices))
    if typeof(C_data) <: Number
        C_data = [C_data]
    end

    # Create and add the new node to the tensor network.
    C = Node(remaining_indices, C_data)
    network.nodes[C_label] = C

    # remove contracted edges
    for index in common_indices
        delete!(network.edges, index)
    end

    # replumb existing edges to point to new node
    # TODO: might become time consuming for large networks
    for index in remaining_indices
        e = network.edges[index]
        if e.src in [A_label, B_label]
            e.src = C_label
        elseif e.dst in [A_label, B_label]
            e.dst = C_label
        end
    end

    # delete the nodes that were contracted
    delete!(network.nodes, A_label)
    delete!(network.nodes, B_label)
    C_label
end




"""
    function compress_tensor_chain!(tng::TensorNetworkCircuit,
                                    nodes::Symbol;
                                    threshold::AbstractFloat=1e-15)

Compress a chain of tensors
"""
function compress_tensor_chain!(tng::TensorNetworkCircuit,
                                nodes::Array{Symbol, 1},
                                executer::AbstractExecuter;
                                threshold::AbstractFloat=1e-15)

    # forward pass
    for i in 1:(length(nodes) - 1)
        left_node = tng.nodes[nodes[i]]
        right_node = tng.nodes[nodes[i+1]]

        left_indices = setdiff(left_node.indices, right_node.indices)
        right_indices = setdiff(right_node.indices, left_node.indices)

        combined_node = contract_pair!(tng, nodes[i], nodes[i+1], executer)

        decompose_tensor!(tng, combined_node, left_indices, right_indices,
                          left_label=nodes[i],
                          right_label=nodes[i+1]
                         )
    end

    # now do a backward pass
    for i in (length(nodes) - 1):-1:1
        left_node = tng.nodes[nodes[i]]
        right_node = tng.nodes[nodes[i+1]]

        left_indices = setdiff(left_node.indices, right_node.indices)
        right_indices = setdiff(right_node.indices, left_node.indices)

        combined_node = contract_pair!(tng, nodes[i], nodes[i+1], executer)

        decompose_tensor!(tng, combined_node, left_indices, right_indices,
                          left_label=nodes[i],
                          right_label=nodes[i+1]
                         )
    end

end
