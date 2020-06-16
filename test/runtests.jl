using QuantExQASM
using PicoQuant

using Compat.Test
using TestSetExtensions
using PyCall

function run_on_qiskit(qasm_str)

    qiskit = pyimport("qiskit")
    circuit = qiskit.QuantumCircuit.from_qasm_str(qasm_str)

    for i in 0:circuit.num_qubits-1
        circuit.measure(i, i)
    end
    
    # Use Aer's qasm_simulator
    simulator = qiskit.Aer.get_backend("qasm_simulator")
    
    # Execute the circuit on the qasm simulator
    job = qiskit.execute(circuit, simulator, shots=1000)
    
    # Grab results from the job
    result = job.result()
    
    # Returns counts
    counts = result.get_counts(circuit)
    return counts
end


function get_statevector_using_picoquant(circ; big_endian=false)
    PicoQuant.InteractiveBackend()
    tn = PicoQuant.convert_qiskit_circ_to_network(circ, decompose=true, transpile=false)
    qubits = circ.n_qubits
    add_input!(tn, "0"^qubits)
    # contract all input nodes together
    full_wf = tn.edges[tn.input_qubits[1]].src
    for i in 2:tn.number_qubits
        full_wf = contract_pair!(tn, full_wf, tn.edges[tn.input_qubits[i]].src)
    end
    plan = PicoQuant.inorder_contraction_plan(tn)
    PicoQuant.contract_network!(tn, plan)
    node = iterate(values(tn.nodes))[1]
    idx_order = [findfirst(x -> x == i, node.indices) for i in tn.output_qubits]
    if !big_endian
        idx_order = idx_order[end:-1:1]
    end
    reshape(permutedims(backend.tensors[:result], idx_order), 2^qubits)
end



function run_on_picoquant(circuit)
    PicoQuant.InteractiveBackend()
    qiskit = pyimport("qiskit")

    cct_s = QuantExQASM.Circuit.to_qasm(circuit, true)

    q_circuit = qiskit.QuantumCircuit.from_qasm_str( String(cct_s) )

    tng = PicoQuant.convert_qiskit_circ_to_network( q_circuit, decompose=false, transpile=false)
    big_endian = true

    qubits = circuit.num_qubits
    add_input!(tng, "0"^qubits)
    plan = inorder_contraction_plan(tng)
    contract_network!(tng, plan)
    node = iterate(values(tng.nodes))[1]
    idx_order = [findfirst(x -> x == i, node.indices) for i in tng.output_qubits]
    if !big_endian
        idx_order = idx_order[end:-1:1]
    end
    reshape(permutedims(PicoQuant.backend.tensors[:result], idx_order), 2^qubits)
end

@testset "All the tests" begin
    @includetests ARGS
end

