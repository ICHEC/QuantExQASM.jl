@testset "Test default NCU operation" begin
    # 4 qubit test
    num_qubits = 4
    cct_s = QuantExQASM.Algorithms.create_ncu_circuit(num_qubits)

    @test begin
        result = run_on_picoquant( cct_s )
        println(result)
        result[16] ≈ 1.0 + 0im
    end
end

@testset "Test NCU+auxiliary (5+3+1) operation" begin
    # 5+3+1 qubit test
    num_qubits = 9
    cct_s = QuantExQASM.Algorithms.create_ncu_circuit(num_qubits, true)

    @test begin
        result = run_on_picoquant( cct_s )
        result[2^Int(ceil(num_qubits/2)) + 2^(num_qubits-1)] ≈ 1.0 + 0im
    end

end