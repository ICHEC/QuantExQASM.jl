@testset "Test default Grover search operation 5 qubits" begin
    # 4 qubit test
    num_qubits = 5
    use_aux_qubits = false
    bit_pattern = 10

    cct = QuantExQASM.Circuit.Circ(num_qubits)
    # Initialise intermediate gates for use in NCU
    QuantExQASM.NCU.init_intermed_gates(cct, num_qubits-1)
    # Run Grover algorithm for given oracle bit-pattern
    QuantExQASM.Grover.run_grover!(cct, collect(0:num_qubits-1), bit_pattern)

    bit_pattern_s = reverse(join(digits(bit_pattern, base=2, pad=num_qubits)))

    # Convert circuit to QASM populated buffer
    cct_s = String(QuantExQASM.Circuit.to_qasm(cct, true) )
    qq_cct = PicoQuant.load_qasm_as_circuit( cct_s )

    @test begin
        counts = run_on_qiskit( cct_s )
        bit_pattern_s in keys(counts)
        counts[bit_pattern_s] > 800
    end

    @test begin
        sv = get_statevector_using_picoquant( qq_cct, big_endian=true )
        rho = abs.(sv).^2
        rho[bit_pattern+1] > 0.8
    end
end

@testset "Test optimised Grover search operation 5+1+3 qubits (Ctrl, Tgt, Aux)" begin
    # 4 qubit test
    num_qubits = 9
    ctrl = 1:5
    tgt = 6
    aux = 7:9
    use_aux_qubits = true
    bit_pattern = 10

    cct = QuantExQASM.Circuit.Circ(num_qubits)
    # Initialise intermediate gates for use in NCU
    QuantExQASM.NCU.init_intermed_gates(cct, num_qubits-1)
    # Run Grover algorithm for given oracle bit-pattern
    QuantExQASM.Grover.run_grover!(cct, collect(vcat(ctrl, tgt)) .- 1, collect(aux) .- 1, bit_pattern)
    # Convert circuit to QASM populated buffer
    cct_s = String(QuantExQASM.Circuit.to_qasm(cct, true))

    bit_pattern_s = reverse(join(digits(bit_pattern, base=2, pad=num_qubits)))

    # Convert circuit to Qiskit circuit
    qq_cct = PicoQuant.load_qasm_as_circuit( cct_s )

    @test begin
        sv = get_statevector_using_picoquant( qq_cct, big_endian=true )
        rho = abs.(sv).^2
        rho[bit_pattern+1] > 0.8
    end
end

@testset "Test default Grover search operation 6 qubits" begin
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
    cct_s = String( QuantExQASM.Circuit.to_qasm(cct, true) )

    bit_pattern_s = reverse(join(digits(bit_pattern, base=2, pad=num_qubits)))

    qq_cct = PicoQuant.load_qasm_as_circuit( cct_s )

    @test begin
        counts = run_on_qiskit( cct_s )

        bit_pattern_s in keys(counts)
        counts[bit_pattern_s] > 800
    end
    @test begin
        sv = get_statevector_using_picoquant( qq_cct, big_endian=true )
        rho = abs.(sv).^2

        rho[bit_pattern+1] > 0.8
    end
end
@testset "Test optimised Grover search operation 6+1+4 qubits (Ctrl, Tgt, Aux)" begin
    # 4 qubit test
    num_qubits = 11
    ctrl = 1:6
    tgt = 7
    aux = 8:11
    use_aux_qubits = true
    bit_pattern = 10

    cct = QuantExQASM.Circuit.Circ(num_qubits)
    # Initialise intermediate gates for use in NCU
    QuantExQASM.NCU.init_intermed_gates(cct, num_qubits-1)
    # Run Grover algorithm for given oracle bit-pattern
    QuantExQASM.Grover.run_grover!(cct, collect(vcat(ctrl, tgt)) .- 1, collect(aux) .- 1, bit_pattern)
    # Convert circuit to QASM populated buffer
    cct_s = String(QuantExQASM.Circuit.to_qasm(cct, true))

    bit_pattern_s = reverse(join(digits(bit_pattern, base=2, pad=num_qubits)))

    # Convert circuit to Qiskit circuit
    qq_cct = PicoQuant.load_qasm_as_circuit( cct_s )

    @test begin
        sv = get_statevector_using_picoquant( qq_cct, big_endian=true )
        rho = abs.(sv).^2
        rho[bit_pattern+1] > 0.8
    end
end