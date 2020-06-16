@testset "Test default NCU operation" begin
    # 4 qubit test
    num_qubits = 4
    use_aux_qubits = false

    # Create empty circuit with qiven qubit count
    cct = QuantExQASM.Circuit.Circ(num_qubits)

    # Initialise default intermediate gates for use in NCU
    QuantExQASM.NCU.init_intermed_gates(cct, num_qubits-1)

    gl = QuantExQASM.GateOps.GateLabel(:x)
    U = QuantExQASM.Circuit.gate_cache[gl]

    ctrl = collect(0:num_qubits-2)
    aux = Int[]
    tgt = num_qubits-1

    for i in ctrl
        QuantExQASM.Circuit.add_gatecall!(cct, QuantExQASM.GateOps.pauli_x(i) )
    end

    QuantExQASM.NCU.apply_ncu!(cct, ctrl, aux, tgt, gl);

    # Convert circuit to QASM populated buffer
    QuantExQASM.Circuit.to_qasm(cct, true, "out.qasm")
    cct_s = QuantExQASM.Circuit.to_qasm(cct, true)

    @test begin
        counts = run_on_qiskit( String(cct_s) )
        "1111" in keys(counts)
        counts["1111"] == 1000
    end
    #@test begin
    #    result = run_on_picoquant( cct )
    #    println(result)
    #    result[16] ≈ 1.0 + 0im
    #end
end

@testset "Test NCU+auxiliary (5+3+1) operation" begin
    # 5+3+1 qubit test
    num_qubits = 9
    use_aux_qubits = true

    ctrl = collect(0:4)
    aux = collect(5:7)
    tgt = 8

    # Create empty circuit with qiven qubit count
    cct = QuantExQASM.Circuit.Circ(num_qubits)

    # Initialise default intermediate gates for use in NCU
    QuantExQASM.NCU.init_intermed_gates(cct, length(ctrl))

    gl = QuantExQASM.GateOps.GateLabel(:x)
    U = QuantExQASM.Circuit.gate_cache[gl]

    for i in ctrl
        QuantExQASM.Circuit.add_gatecall!(cct, QuantExQASM.GateOps.pauli_x(i) )
    end

    @test begin
        result = run_on_picoquant( cct )
        result[2^length(ctrl)] ≈ 1.0 + 0im
    end

    QuantExQASM.NCU.apply_ncu!(cct, ctrl, aux, tgt, gl);

    # Convert circuit to QASM populated buffer
    QuantExQASM.Circuit.to_qasm(cct, true, "out.qasm")
    cct_s = QuantExQASM.Circuit.to_qasm(cct, true)

    @test begin
        result = run_on_picoquant( cct )
        result[2^length(ctrl) + 2^tgt] ≈ 1.0 + 0im
    end

end