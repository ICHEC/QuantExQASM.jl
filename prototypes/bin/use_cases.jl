using Pkg
using Logging

# Activate local packages from repo if not installed
try
    using QuantExQASM
catch
    mod_path = normpath(joinpath(@__DIR__, "..", "QuantExQASM"))
    Pkg.activate(mod_path) 
    using QuantExQASM
end
try
    using PicoQuant
catch
    mod_path = normpath(joinpath(@__DIR__, "..", "PicoQuant"))
    Pkg.activate(mod_path) 
    using PicoQuant
end
using ArgParse
using PyCall

# *************************************************************************** #
#                 Use-case circuit generation functions
# *************************************************************************** #

"""
    function create_grover_circuit(qubits::Integer)
"""
function create_grover_circuit(num_qubits::Int)

    qubit_indices = collect(0:num_qubits-1)
    cl_indices = qubit_indices # Measure all

    q_reg_label = "qu"
    c_reg_label = "cl"

    marked_state = 5

    # Create Grover operation 
    cct = QuantExQASM.GateOps.gen_qasm_header(num_qubits, length(cl_indices), q_reg_label, c_reg_label)
    cct *= QuantExQASM.Grover.run_grover(q_reg_label, qubit_indices, marked_state)
    cct *= QuantExQASM.GateOps.measure_qubits_all(qubit_indices, cl_indices, q_reg_label, c_reg_label)
    return PicoQuant.load_qasm_as_circuit(cct)
end

# *************************************************************************** #

"""
    function create_phase_oracle_circuit(qubits::Integer)
"""
function create_phase_oracle_circuit(num_qubits::Int)

    qubit_indices = collect(0:num_qubits-1)
    cl_indices = qubit_indices # Measure all

    q_reg_label = "qu"
    c_reg_label = "cl"

    marked_state = 1

    # Create Grover operation 
    cct = QuantExQASM.GateOps.gen_qasm_header(num_qubits, length(cl_indices), q_reg_label, c_reg_label)
    cct *= QuantExQASM.Oracle.bitstring_phase_oracle(q_reg_label, marked_state, qubit_indices[1:end-1], qubit_indices[end])
    cct *= QuantExQASM.GateOps.measure_qubits_all(qubit_indices, cl_indices, q_reg_label, c_reg_label)
    
    return PicoQuant.load_qasm_as_circuit(cct)
end

# *************************************************************************** #

"""
    function create_ncu_circuit(qubits::Integer)

Generate n-controlled unitary circuit acting on given number of qubits.

Defaults to |11...110> register, and applies nCX to generate |11..111>
"""
function create_ncu_circuit(num_qubits::Int)

    X = QuantExQASM.NCU.GateOps.default_gates["X"]

    use_aux_qubits = true

    if use_aux_qubits == true
        ctrl = collect(0:num_qubits-2)
        tgt = num_qubits-1
        aux = collect( num_qubits: num_qubits + length(ctrl) - 3)
    else
        ctrl = collect(0:num_qubits-2)
        aux = Array
        tgt = length(ctrl)
    end

    num_classical = num_qubits # Measure only the target line
    q_reg_label = "qu"
    c_reg_label = "cl"

    # Create NCU operation 
    cct = QuantExQASM.GateOps.gen_qasm_header(2*num_qubits-3, num_classical, q_reg_label, c_reg_label)
    for i in ctrl
        cct *= QuantExQASM.GateOps.apply_gate_x(q_reg_label, i)
    end
    cct *= QuantExQASM.NCU.apply_ncu(q_reg_label, ctrl, aux, tgt, X, 0)
    cct *= QuantExQASM.GateOps.measure_qubits_all(vcat(ctrl,tgt), vcat(ctrl,tgt), q_reg_label, c_reg_label)
    return PicoQuant.load_qasm_as_circuit(cct)
end

# *************************************************************************** #

"""
    function create_qft_circuit(qubits::Integer)

Generate QFT circuit acting on given number of qubits
"""
function create_qft_circuit(qubits::Integer, qiskit_generated::Bool=false)
    qiskit = pyimport("qiskit")
    qfts = pyimport("qiskit.aqua.components.qfts")

    if !qiskit_generated
        circ = qiskit.QuantumCircuit(qubits)
        for i in 1:qubits
            circ.h(i-1)
            for j in i:qubits-1
                circ.cu1(π/(2^(j-i+2)), j, i-1)
            end
            circ.barrier()
        end
        for i = 1:convert(Integer, floor(qubits//2))
            circ.swap(i-1, qubits -i)
        end
    else
        circ = qiskit.QuantumCircuit(qubits)
        circ = qfts.Standard(qubits).construct_circuit(mode="circuit",
                                                       qubits=circ.qregs[1],
                                                       circuit=circ)
    end

    return circ
end

# *************************************************************************** #
#                           I/O Functions
# *************************************************************************** #

function latex_output(filename, circuit)
    latex_filename = "$filename.tex"
    circuit.draw(output="latex_source", filename=latex_filename)
end
function text_output(filename, circuit)
    text_filename = "$filename.txt"
    circuit.draw(output="text", filename=text_filename)
end
function tn_output(filename, circuit)
    tng = PicoQuant.convert_to_tensor_network_graph(circuit)
    tn_filename = "$filename.json"
    open(tn_filename, "w") do io
        write(io, PicoQuant.to_json(tng, parsed_args["indent"]))
    end
end
function qasm_output(filename, circuit)
    qasm_filename = "$filename.qasm"
    open(qasm_filename, "w") do io
        write(io, circuit.qasm())
    end
end

output_fcall = Dict(
    "qasm"=>qasm_output,
    "tn"=>tn_output,
    "txt"=>text_output,
    "latex"=>latex_output,
)

# *************************************************************************** #
#                       Parsing and main functions
# *************************************************************************** #

"""
    function parse_commandline()

Parse command line options and return argument dictionary
"""
function parse_commandline()
    s = ArgParseSettings()
        @add_arg_table! s begin
            "--verbose", "-v"
                help = "Enable verbose logging output."
                action = :store_true
            "--use-case", "-u"
                help = "Use-case to generate: qft, ncu, grover"
                arg_type = String
                default = "qft"
            "--gen-type", "-g"
                help = "Output format(s) to generate for. Options: \"qasm,latex,tn,txt\". Use comma between multiple options."
                arg_type = String
                default = "qasm"
            "--number-qubits", "-n"
                help = "Size of circuit"
                arg_type = Int
                required = true
            "--output", "-o"
                help = "Output file label. Format will be <output>.<gen-type>"
                arg_type = String
                default = ""
            "--latex-circuit"
                help = "Generate latex for circuit diagram"
                action = :store_true
            "--text-circuit"
                help = "Generate ASCII for circuit diagram"
                action = :store_true
            "--qiskit-generated"
                help = "Use qiskit's inbuilt QFT class instead of constructing circuit manually "
                action = :store_true
        end
        return parse_args(s)
end

# *************************************************************************** #

function main()
    parsed_args = parse_commandline()
    # ************************* #
    verbosity = parsed_args["verbose"]
    if verbosity == true
        io = open("out.log", "w+")
        logger = SimpleLogger(io, Logging.Debug)
        global_logger(logger)
    end

    use_case = parsed_args["use-case"]
    qubits = parsed_args["number-qubits"]
    # ************************* #
    # ************************* #
    if use_case == "qft"
        circuit= create_qft_circuit(qubits, parsed_args["qiskit-generated"])
    elseif use_case == "ncu"
        circuit = create_ncu_circuit(qubits)
    elseif use_case == "grover"
        circuit = create_grover_circuit(qubits)
    elseif use_case == "oracle"
        circuit = create_phase_oracle_circuit(qubits)
    else
        throw("Unknown use-case")
    end
    # ************************* #
    # ************************* #
    if parsed_args["output"] == ""
        filename = "$(use_case)_$(qubits)"
    else
        filename = parsed_args["output"]
    end
    # ************************* #
    # ************************* #

    out_types = split(parsed_args["gen-type"],",")
    for i in out_types
        output_fcall[i](filename, circuit)
    end
end

# *************************************************************************** #
# *************************************************************************** #

main()
