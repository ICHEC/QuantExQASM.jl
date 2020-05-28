using Pkg
using Logging
using UnicodePlots

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
    function create_grover_circuit(qubits::Integer, bit_pattern::Int=11)
"""
function create_grover_circuit(num_qubits::Int, bit_pattern::Int=11)
    # Create empty circuit with qiven qubit count
    cct = QuantExQASM.Circuit.Circ(num_qubits)
    # Initialise intermediate gates for use in NCU
    QuantExQASM.NCU.init_intermed_gates(cct, num_qubits-1)
    # Run Grover algorithm for given oracle bit-pattern
    QuantExQASM.Grover.run_grover!(cct, collect(0:num_qubits-1), bit_pattern)
    # Convert circuit to QASM populated buffer
    cct_s = QuantExQASM.Circuit.to_qasm(cct, true)
    # Convert to Qiskit loaded qasm circuit
    #println( String(cct_s) )
    #exit()
    return PicoQuant.load_qasm_as_circuit( String(cct_s) )
end


# *************************************************************************** #

"""
    function create_phase_oracle_circuit(qubits::Integer, bit_pattern::Int=11)
"""
function create_phase_oracle_circuit(num_qubits::Int, bit_pattern::Int=11)
    # Create empty circuit with qiven qubit count
    cct = QuantExQASM.Circuit.Circ(num_qubits)
    # Initialise intermediate gates for use in NCU
    QuantExQASM.NCU.init_intermed_gates(cct, num_qubits-1)
    # Run Grover algorithm for given oracle bit-pattern
    QuantExQASM.Oracle.bitstring_phase_oracle!(cct, bit_pattern, collect(0:num_qubits-2), num_qubits-1)
    # Convert circuit to QASM populated buffer
    cct_s = QuantExQASM.Circuit.to_qasm(cct, true)
    # Convert to Qiskit loaded qasm circuit
    return PicoQuant.load_qasm_as_circuit( String(cct_s) )
end
# *************************************************************************** #

"""
    function create_ncu_circuit(qubits::Integer)

Generate n-controlled unitary circuit acting on given number of qubits.

Defaults to |11...110> register, and applies nCX to generate |11..111>
"""

function create_ncu_circuit(num_qubits::Int, gl::QuantExQASM.GateOps.GateLabel=QuantExQASM.GateOps.GateLabel(:x), use_aux_qubits::Bool=true)

    # Create empty circuit with qiven qubit count
    cct = QuantExQASM.Circuit.Circ(num_qubits)

    # Initialise default intermediate gates for use in NCU
    QuantExQASM.NCU.init_intermed_gates(cct, num_qubits-1)
    
    U = QuantExQASM.Circuit.gate_cache[gl]

    if use_aux_qubits == true && (num_qubits/2+1 >= 4)
        ctrl = collect(0:Int( floor((num_qubits-1)//2) ))
        aux = collect(1 + Int( floor((num_qubits-1)//2)):num_qubits-2)
        tgt = num_qubits-1
    else
        ctrl = collect(0:num_qubits-2)
        aux = Int[]
        tgt = num_qubits-1
    end

    for i in ctrl
        QuantExQASM.Circuit.add_gatecall!(cct, QuantExQASM.GateOps.pauli_x(i) )
    end

    QuantExQASM.NCU.apply_ncu!(cct, ctrl, aux, tgt, gl);

    # Convert circuit to QASM populated buffer
    QuantExQASM.Circuit.to_qasm(cct, true, "out.qasm")
    cct_s = QuantExQASM.Circuit.to_qasm(cct, true)

    return PicoQuant.load_qasm_as_circuit( String(cct_s) )
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
    tng = PicoQuant.convert_qiskit_circ_to_network(circuit)
    tn_filename = "$filename.json"
    open(tn_filename, "w") do io
        write(io, PicoQuant.to_json(tng, 4))
    end
end
function qasm_output(filename, circuit)
    qasm_filename = "$filename.qasm"

    open(qasm_filename, "w") do io
        write(io, circuit.decompose().qasm())
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
            "--run-qiskit"
                help = "Run the example circuit on Qiskit's QASM simulator"
                action = :store_true
        end
        return parse_args(s)
end

# *************************************************************************** #

function run_on_qiskit(circuit)

    qiskit = pyimport("qiskit")

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

function run_on_qiskit2(circuit)
    qiskit = pyimport("qiskit")

    for i in 0:circuit.num_qubits-1
        circuit.measure(i, i)
    end

    backend = qiskit.providers.aer.QasmSimulator()
    backend_options = Dict("method"=>"matrix_product_state", "max_parallel_threads"=>2, "max_parallel_shots"=>4)
    
    # Circuit execution
    job = qiskit.execute(circuit, backend, shots=1000, backend_options=backend_options)
    # Grab results from the job
    result = job.result()
    
    # Returns counts
    counts = result.get_counts(circuit)
    return counts 
end

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
    if parsed_args["run-qiskit"] == true
        counts = run_on_qiskit(circuit)
        println(counts)

        k = collect(keys(counts))
        v = convert(Array{Int}, collect(values(counts)))
        b = barplot(k, v)
        print(b)
        return;
    end
    # ************************* #
    if parsed_args["output"] == ""
        filename = "$(use_case)_$(qubits)"
    else
        filename = parsed_args["output"]
    end
    # ************************* #
    # ************************* #
    #Output multiple types by giving csv labels
    out_types = split(parsed_args["gen-type"],",")
    for i in out_types
        output_fcall[i](filename, circuit)
    end
end

# *************************************************************************** #
# *************************************************************************** #

main()
