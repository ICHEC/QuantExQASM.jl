@testset "Test NCU operation" begin
    # 4 qubit test
    num_qubits = 4
    use_aux_qubits = false

    # Create empty circuit with qiven qubit count
    cct = QuantExQASM.Circuit.Circ(num_qubits)

    # Initialise default intermediate gates for use in NCU
    QuantExQASM.NCU.init_intermed_gates(cct, num_qubits-1)

    gl = QuantExQASM.GateOps.GateLabel(:x)

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
    @test begin
        counts = run_on_qiskit( String(cct_s) )
        "1111" in keys(counts)
        counts["1111"] == 1000
    end
    @test begin
        result = run_on_picoquant( cct )
        println(result)
        result[16] ≈ 1.0 + 0im
    end

end

