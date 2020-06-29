module Algorithms
export create_ncu_circuit, create_grover_circuit

include("algs/NCU.jl")
include("algs/Oracle.jl")
include("algs/Diffusion.jl")
include("algs/Grover.jl")

using QuantExQASM.GateOps
using QuantExQASM.Circuit
using PicoQuant

function create_ncu_circuit(num_qubits::Int, use_aux::Bool=false, gate_label::GateOps.GateLabel=GateOps.GateLabel(:x))
    # Create empty circuit with qiven qubit count
    cct = Circuit.Circ(num_qubits)

    # Initialise default intermediate gates for use in NCU
    NCU.init_intermed_gates(cct, num_qubits-1)

    try # check if accessing Gate exists
        U = Circuit.gate_cache[gate_label]
    catch e
        println("GateLabel $(gate_label) does not exist.")
        throw(e)
    end

    if use_aux == false || num_qubits < 9 # If insufficient qubits, use unoptimised NCU
        ctrl = collect(0:num_qubits-2)
        aux = Int[]
        tgt = num_qubits-1
    else
        num_ctrl = Int(ceil(num_qubits/2))
        num_aux = num_qubits - num_ctrl -1
        ctrl = collect(0:num_ctrl-1)
        aux = collect(num_ctrl:num_ctrl+num_aux-1)
        tgt = num_qubits-1
    end

    for i in ctrl
        Circuit.add_gatecall!(cct, GateOps.pauli_x(i) )
    end

    NCU.apply_ncu!(cct, ctrl, aux, tgt, gate_label);

    # Convert circuit to QASM populated buffer
    Circuit.to_qasm(cct, true, "out.qasm")
    cct_s = String( Circuit.to_qasm(cct, true) )
    qq_cct = PicoQuant.load_qasm_as_circuit( cct_s )
    return qq_cct
end

function create_grover_circuit(num_qubits::Int, use_aux::Bool=false, bit_pattern::Int=11)
    # Create empty circuit with qiven qubit count
    cct = Circuit.Circ(num_qubits)

    # Initialise default intermediate gates for use in NCU
    NCU.init_intermed_gates(cct, num_qubits-1)

    if use_aux == false || num_qubits < 9 # If insufficient qubits, use unoptimised NCU
        Grover.run_grover!(cct, collect(0:num_qubits-1), bit_pattern)
    else
        num_ctrl = Int(ceil(num_qubits/2))
        num_aux = num_qubits - num_ctrl -1
        ctrl = collect(0:num_ctrl-1)
        aux = collect(num_ctrl:num_ctrl+num_aux-1)
        tgt = num_qubits-1

        Grover.run_grover!(cct, vcat(ctrl, tgt), aux, bit_pattern)
    end

    #bit_pattern_s = reverse(join(digits(bit_pattern, base=2, pad=num_qubits)))

    cct_s = String(Circuit.to_qasm(cct, true) )
    qq_cct = PicoQuant.load_qasm_as_circuit( cct_s )
    return qq_cct
end

end