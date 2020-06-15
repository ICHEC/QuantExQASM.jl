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

    bit_pattern_s = reverse(join(digits(bit_pattern, base=2, pad=num_qubits)))

    # Convert circuit to QASM populated buffer
    cct_s = QuantExQASM.Circuit.to_qasm(cct, true)

    @test begin
        counts = run_on_qiskit( String(cct_s) )
        bit_pattern_s in keys(counts)
        counts[bit_pattern_s] > 800
    end

    @test begin
        sv = run_on_picoquant( cct )
        rho = abs.(sv).^2
        rho[bit_pattern+1] > 0.8
    end
end

@testset "Test Grover search operation 6 qubits" begin
    # 7 qubit test
    num_qubits = 6
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

    @test begin
        counts = run_on_qiskit( String(cct_s) )

        bit_pattern_s in keys(counts)
        counts[bit_pattern_s] > 800
    end
    @test begin
        sv = run_on_picoquant( cct )
        rho = abs.(sv).^2
        rho[bit_pattern+1] > 0.8
    end
end


