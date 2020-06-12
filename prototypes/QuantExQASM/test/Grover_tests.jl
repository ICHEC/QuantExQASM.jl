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

@testset "Test Grover search operation 5 qubits" begin
    # 4 qubit test
    num_qubits = 5
    use_aux_qubits = false
    bit_pattern = 10

    cct = QuantExQASM.Circuit.Circ(num_qubits)
    # Initialise intermediate gates for use in NCU
    QuantExQASM.NCU.init_intermed_gates(cct, num_qubits-1)
    # Run Grover algorithm for given oracle bit-pattern
    QuantExQASM.Grover.run_grover!(cct, collect(0:num_qubits-1), bit_pattern)
    # Convert circuit to QASM populated buffer
    cct_s = QuantExQASM.Circuit.to_qasm(cct, true)

    bit_pattern_s = reverse(join(digits(bit_pattern, base=2,pad=num_qubits)))

    # Convert circuit to QASM populated buffer
    cct_s = QuantExQASM.Circuit.to_qasm(cct, true)

    counts = run_on_qiskit( String(cct_s) )

    @test begin
        bit_pattern_s in keys(counts)
        counts[bit_pattern_s] > 800
    end
end


@testset "Test Grover search operation 7 qubits" begin
    # 7 qubit test
    num_qubits = 7
    use_aux_qubits = false
    bit_pattern = 37

    cct = QuantExQASM.Circuit.Circ(num_qubits)
    # Initialise intermediate gates for use in NCU
    QuantExQASM.NCU.init_intermed_gates(cct, num_qubits-1)
    # Run Grover algorithm for given oracle bit-pattern
    QuantExQASM.Grover.run_grover!(cct, collect(0:num_qubits-1), bit_pattern)
    # Convert circuit to QASM populated buffer
    cct_s = QuantExQASM.Circuit.to_qasm(cct, true)

    bit_pattern_s = reverse(join(digits(bit_pattern, base=2, pad=num_qubits)))

    # Convert circuit to QASM populated buffer
    cct_s = QuantExQASM.Circuit.to_qasm(cct, true)

    counts = run_on_qiskit( String(cct_s) )

    @test begin
        bit_pattern_s in keys(counts)
        counts[bit_pattern_s] > 800
    end
end


