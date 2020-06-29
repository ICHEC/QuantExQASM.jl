@testset "Test default NCU operation" begin
    # 5 qubit test
    num_qubits = 5
    cct_s = QuantExQASM.Algorithms.create_ncu_circuit(num_qubits)

    @test begin
        result = get_statevector_using_picoquant( cct_s, big_endian=true )
        result[2^num_qubits] ≈ 1.0 + 0im
    end
    @test begin
        result = get_statevector_using_picoquant_mps( cct_s )
        result[2^num_qubits] ≈ 1.0 + 0im
    end
end

@testset "Test NCU+auxiliary (5+3+1) operation" begin
    # 5+3+1 qubit test
    num_qubits = 9
    cct_s = QuantExQASM.Algorithms.create_ncu_circuit(num_qubits, true)

    @test begin
        result = get_statevector_using_picoquant( cct_s, big_endian=true )
        result[2^Int(ceil(num_qubits/2)) + 2^(num_qubits-1)] ≈ 1.0 + 0im
    end
    @test begin
        result = get_statevector_using_picoquant_mps( cct_s )
        result[2^Int(ceil(num_qubits/2)) + 2^(num_qubits-1)] ≈ 1.0 + 0im
    end
end