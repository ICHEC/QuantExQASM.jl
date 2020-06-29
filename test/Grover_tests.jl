@testset "Test default Grover search operation 5 qubits" begin
    num_qubits = 5
    use_aux_qubits = false
    bit_pattern = 10

    qq_cct = QuantExQASM.Algorithms.create_grover_circuit(num_qubits, use_aux_qubits, bit_pattern)

    @test begin
        sv = get_statevector_using_picoquant( qq_cct, big_endian=true )
        rho = abs.(sv).^2
        rho[bit_pattern+1] > 0.8
    end
end

@testset "Test optimised Grover search operation 5+1+3 qubits (Ctrl, Tgt, Aux)" begin
    num_qubits = 9
    use_aux_qubits = true
    bit_pattern = 10

    # Convert circuit to Qiskit circuit
    qq_cct = QuantExQASM.Algorithms.create_grover_circuit(num_qubits, use_aux_qubits, bit_pattern)

    @test begin
        sv = get_statevector_using_picoquant( qq_cct, big_endian=true )
        rho = abs.(sv).^2
        rho[bit_pattern+1] > 0.8
    end
end

@testset "Test default Grover search operation 6 qubits" begin
    num_qubits = 6
    use_aux_qubits = false
    bit_pattern = 37

    qq_cct = QuantExQASM.Algorithms.create_grover_circuit(num_qubits, use_aux_qubits, bit_pattern)

    @test begin
        sv = get_statevector_using_picoquant( qq_cct, big_endian=true )
        rho = abs.(sv).^2
        rho[bit_pattern+1] > 0.8
    end
end

@testset "Test optimised Grover search operation 6+1+4 qubits (Ctrl, Tgt, Aux)" begin
    num_qubits = 11
    use_aux_qubits = true
    bit_pattern = 10

    qq_cct = QuantExQASM.Algorithms.create_grover_circuit(num_qubits, use_aux_qubits, bit_pattern)

    @test begin
        sv = get_statevector_using_picoquant( qq_cct, big_endian=true )
        rho = abs.(sv).^2
        rho[bit_pattern+1] > 0.8
    end
end