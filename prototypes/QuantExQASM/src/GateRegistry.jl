module GateRegistry

import GateOps
GateLabel = GateOps.GateLabel
Gate = GateOps.Gate

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#                        Map GateLabel to Gates
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

export pauli_x, pauli_y, pauli_z, hadamard, u
export r_x, r_y, r_z, r_phase, swap
export c_pauli_x, c_pauli_y, c_pauli_z, c_u
export c_r_x, c_r_y, c_r_z, c_r_phase

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

# For each unique GateLabel we require an assigned matrix value
GateMap = Dict{GateLabel, Gate}(
    GateLabel(:x)=>Gate([0 1; 1 0] .+ 0im),
    GateLabel(:y)=>Gate([0 -1im; 1im 0]),
    GateLabel(:z)=>Gate([1 0; 0 -1] .+ 0im),
    GateLabel(:h)=>Gate([1 1; 1 -1]*(1/sqrt(2)) .+ 0im)
)

# List available operators
# Generalise later to accept various types
Operators = Dict{String, Function}(
    "x"=>GateOps.pauli_x,
    "y"=>GateOps.pauli_y,
    "z"=>GateOps.pauli_z,
    "h"=>GateOps.hadamard,
    "u"=>GateOps.u,
    "r_x"=>GateOps.r_x,
    "r_y"=>GateOps.r_y,
    "r_z"=>GateOps.r_z,
    "r_phase"=>GateOps.r_phase,
    "c_x"=>GateOps.c_pauli_x,
    "c_y"=>GateOps.c_pauli_y,
    "c_z"=>GateOps.c_pauli_z,
    "c_u"=>GateOps.c_u,
    "c_r_x"=>GateOps.c_r_x,
    "c_r_y"=>GateOps.c_r_y,
    "c_r_z"=>GateOps.c_r_z,
    "c_r_phase"=>GateOps.c_r_phase
)

# List available modifications to gate
GateModifiers = Dict{String, Function}(
    "sqrt"=>GateOps.sqrt,
    "adj"=>GateOps.adjoint,
    "e2m"=>GateOps.euler_to_mat,
    "m2e"=>GateOps.mat_to_euler,
    "no_op"=>x->x
)

struct Registry
    map_label_gate::Dict{GateLabel, Gate}
    operators::Dict{String, Function}
    gate_mods::Dict{String, Function}
end
GateReg = Registry(GateMap, Operators, GateModifiers)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#                               Get and cache
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
@memoize function apply_gate(gate_func::Function, gate_modifier::Union{Nothing, Function}, params::Tuple)
    if gate_modifier == nothing
        return gate_func(params...)
    else
        gate_modifier(gate_func(params...)) # Construct new gate, symbol, label, and cache before returning
    end
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#                               Gate cache
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

function create_reg_entity(label::String, mat::Matrix{<:Number})
    g = Gate(mat)
    l = GateLabel(label)
    return g=>l
end

function create_reg_operator(label::String, f::Function, params::Tuple)
    l = GateLabel(label)
    g = f(params...)
    add_gate!(l, g)
    push!(Operators, )
    return g=>l
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#                          Boilerplate functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
function add_gate!(g_map::Dict{GateLabel, Gate}, label::GateLabel, gate::Gate)
    if haskey(g_map, label)
        return
    else
        push!(g_map, label=>gate)
    end
end
function add_gate!(label::GateLabel, gate::Gate)
    if haskey(GateMap, label)
        return
    else
        push!(GateMap, label=>gate)
    end
end
function add_gate!(label_gate::Pair{GateLabel, Gate})
    if haskey(GateMap, label_gate.first)
        return
    else
        push!(GateMap, label_gate)
    end
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

function remove_gate!(g_map::Dict{GateLabel, Gate}, label::GateLabel)
    delete!(g_map, label)
end

function remove_gate!(label::GateLabel)
    delete!(GateMap, label)
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

function get_gate(label::GateLabel)
    if ~haskey(GateMap, label)
        return GateMap[label]
    else
        println("Key $(label) does not exist")
    end
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

function has_gate(label::GateLabel)
    return haskey(GateMap, label)
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

function list_gates(;val=false)
    if ~val
        return keys(GateMap)
    else
        return values(GateMap)
    end
end

function list_gate_map(;val=false)    
    return GateMap
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

@memoize function Base.:sqrt(g::GateOps.GateLabelP)
    m = sqrt(GateMap[g])
    GateOps.GateLabelP
end

function mat_to_euler(gate::GateOps.GateLabel)
    return mat_to_euler(GateMap[gate])
end
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

end