module NCU
using QuantExQASM.GateOps
using QuantExQASM.Circuit

ops_mat_cache = Dict{GateOps.GateLabel, Matrix{<:Number}}()

# =========================================================================== #
# 
# =========================================================================== #

function init_intermed_gates(circ::Circuit.Circ, num_ctrl::Union{Nothing, Int})
    for k in circ.gate_set
        gen_intermed_gates(num_ctrl == nothing ? 8 : num_ctrl, k)
    end
end

function register_gate(circ::Circuit.Circ, U::GateOps.GateLabel, gate::Matrix{<:Number})
    # Some characters are reserved for internal gates
    @assert !(String(U.label) in ["x","y","z","h","s","t","a","c"])
    Circuit.gate_cache[U] = gate
end

function gen_intermed_gates(ctrl_depth::Int, U::GateOps.GateLabel)
    su,asu = get_intermed_gate(U)
    for i in 2:ctrl_depth-1
        su,asu = get_intermed_gate(su)
    end
end

function get_intermed_gate(U::GateOps.GateLabel)
    su = GateOps.GateLabel( Symbol("s" * String(U.label)) )
    asu = GateOps.GateLabel( Symbol("as" * String(U.label)) )

    if haskey(Circuit.gate_cache, su)
        SU = Circuit.gate_cache[su]
        ASU = Circuit.gate_cache[asu]
    else

        SU = sqrt(Circuit.gate_cache[U])
        ASU = collect(adjoint(sqrt(Circuit.gate_cache[U])))

        # Cache the gates
        Circuit.gate_cache[su] = SU
        Circuit.gate_cache[asu] = ASU

        # Print the newly generated gates to qasm file for use.
        # Gates can be defined anywhere in file.

        #GateOps.create_gates_nonparam(su, 1)
        #GateOps.create_gates_nonparam(asu, 1)
        #GateOps.create_gates_nonparam(su, 2)
        #GateOps.create_gates_nonparam(asu, 2)
    end

    return su,asu
end

# =========================================================================== #
# 
# =========================================================================== #

function apply_cx!(c::Circuit.Circ, ctrl, tgt, reg) 
    Circuit.add_gatecall!(c, GateOps.c_pauli_x(ctrl, tgt, reg))
end

function apply_cu!(c::Circuit.Circ, ctrl, tgt, reg, gl::GateOps.GateLabel)
    if String(gl.label)[end] == 'x'
        glz_s = String(gl.label)[1:end-1] * 'z'
        gl_z = GateOps.GateLabel(Symbol(glz_s))

        push!(c.gate_set, gl_z)

        Circuit.add_gatecall!(c, GateOps.hadamard(tgt, reg))
        Circuit.add_gatecall!(c, GateOps.c_u(gl_z, ctrl, tgt, GateOps.u( GateOps.GateLabel(Symbol(split(glz_s, "c_")[end])), tgt, reg), reg))
        Circuit.add_gatecall!(c, GateOps.hadamard(tgt, reg))

    elseif String(gl.label)[end] == 'z'
        Circuit.add_gatecall!(c, GateOps.c_u(gl, ctrl, tgt, GateOps.pauli_z(tgt, reg), reg))
    else
        error("Currently only PauliX and PauliZ decomposed gates supported")
    end
end

# =========================================================================== #
# 
# =========================================================================== #

"""
    apply_ncu(ctrls::Vector{Int}, aux::Vector{Int}, tgt::Int, U::GateOps.GateLabel)

Apply an n-qubit controlled gate operation on the given target qubit.
Ensure the gate corresponding with symbol U is registered with g_map before use.

# Arguments
- `ctrls::Vector{Int}`: 
- `aux::Vector{Int}`: 
- `tgt::Int`:
- `U::Symbol`:
"""
function apply_ncu(q_ctrl::Vector, q_aux::Vector, q_tgt, U::GateOps.GateLabel)
    cct = Circuit.Circ()

    # Check global circuit cache for gate
    if ~haskey(Circuit.gate_cache, U )
        error("Gate $(U) does not exist")
    end

    su,asu = get_intermed_gate(U)

    if length(q_ctrl) == 2
#=
        cct.add_gatecall( GateOps.c_u(su, q_tgt, q_ctrl[2]) )
        cct.add_gatecall( GateOps.c_pauli_x(q_ctrl[2], q_ctrl[1]) )
        cct.add_gatecall( GateOps.c_u(asu, q_tgt, q_ctrl[2]) )
        cct.add_gatecall( GateOps.c_pauli_x(q_ctrl[2], q_ctrl[1]) )
        cct.add_gatecall( GateOps.c_u(su, q_tgt, q_ctrl[1]) )
=#

        apply_cu!(cct, q_ctrl[2], q_tgt, nothing, su)
        apply_cx!(cct, q_ctrl[1], q_ctrl[2], nothing)
        apply_cu!(cct, q_ctrl[2], q_tgt,  nothing, asu)
        apply_cx!(cct, q_ctrl[1], q_ctrl[2], nothing)
        apply_cu!(cct, q_ctrl[1], q_tgt,  nothing, su)

    elseif length(q_ctrl) == 3
        #ssu = sqrt(su)
        #assu = adjoint(ssu)
        ssu,assu = get_intermed_gate(su)
#=
        Circuit.add_gatecall!(cct, GateOps.c_u(ssu, q_tgt, q_ctrl[1]) )
        Circuit.add_gatecall!(cct, GateOps.c_pauli_x(q_ctrl[2], q_ctrl[1]) )
        Circuit.add_gatecall!(cct, GateOps.c_u(assu, q_tgt, q_ctrl[2]) )
        Circuit.add_gatecall!(cct, GateOps.c_pauli_x(q_ctrl[2], q_ctrl[1]) )
        Circuit.add_gatecall!(cct, GateOps.c_u(ssu, q_tgt, q_ctrl[2]) )
        Circuit.add_gatecall!(cct, GateOps.c_pauli_x(q_ctrl[3], q_ctrl[2]) )
        Circuit.add_gatecall!(cct, GateOps.c_u(assu, q_tgt, q_ctrl[3]) )
        Circuit.add_gatecall!(cct, GateOps.c_pauli_x(q_ctrl[3], q_ctrl[1]) )
        Circuit.add_gatecall!(cct, GateOps.c_u(ssu, q_tgt, q_ctrl[3]) )
        Circuit.add_gatecall!(cct, GateOps.c_pauli_x(q_ctrl[3], q_ctrl[2]) )
        Circuit.add_gatecall!(cct, GateOps.c_u(assu, q_tgt, q_ctrl[3]) )
        Circuit.add_gatecall!(cct, GateOps.c_pauli_x(q_ctrl[3], q_ctrl[1]) )
        Circuit.add_gatecall!(cct, GateOps.c_u(ssu, q_tgt, q_ctrl[3]) )
=#
        apply_cu!(cct, q_ctrl[1], q_tgt, nothing, ssu)
        apply_cx!(cct, q_ctrl[1], q_ctrl[2], nothing)
        apply_cu!(cct, q_ctrl[2], q_tgt, nothing, assu)
        apply_cx!(cct, q_ctrl[1], q_ctrl[2], nothing)
        apply_cu!(cct, q_ctrl[2], q_tgt, nothing, ssu)
        apply_cx!(cct, q_ctrl[2], q_ctrl[3], nothing)
        apply_cu!(cct, q_ctrl[3], q_tgt, nothing, assu)
        apply_cx!(cct ,q_ctrl[1], q_ctrl[3], nothing)
        apply_cu!(cct, q_ctrl[3], q_tgt, nothing, ssu)
        apply_cx!(cct, q_ctrl[2], q_ctrl[3], nothing)
        apply_cu!(cct, q_ctrl[3], q_tgt, nothing, assu)
        apply_cx!(cct, q_ctrl[1], q_ctrl[3], nothing)
        apply_cu!(cct, q_ctrl[3], q_tgt, nothing, ssu)
        
    elseif (length(q_ctrl)>=5) && (length(q_aux)>=length(q_ctrl)-2) && (U == GateOps.GateLabel(:x))
        append!(cct, apply_ncu([q_ctrl[end], q_aux[ 1 + length(q_ctrl)-3 ]], Int[], q_tgt, U) )
        
        for i in reverse(2:length(q_ctrl)-2)
            append!(cct, apply_ncu([q_ctrl[1 + i], q_aux[1 + (i-2)]], Int[], q_aux[1 + (i-1)], U) )
        end
        
        append!(cct, apply_ncu([q_ctrl[1], q_ctrl[2]], Int[], q_aux[1], U) )
        for i in 2:length(q_ctrl)-2
            append!(cct, apply_ncu([q_ctrl[1 + i], q_aux[1 + (i-2)]], Int[], q_aux[1 + (i-1)], U) )
        end
        
        append!(cct, apply_ncu([q_ctrl[end], q_aux[1 + length(q_ctrl) - 3]], Int[], q_tgt, U) )
        for i in reverse(2:length(q_ctrl)-2)
            append!(cct, apply_ncu([q_ctrl[1 + i], q_aux[1 + (i-2)]], Int[], q_aux[1 + (i-1)], U) )
        end
        
        append!(cct, apply_ncu([q_ctrl[1], q_ctrl[2]], Int[], q_aux[1], U) )
        for i in 2:length(q_ctrl)-2
            append!(cct, apply_ncu([q_ctrl[1 + i], q_aux[1 + (i-2)]], Int[], q_aux[1 + (i-1)], U) )
        end

    else
#=
        Circuit.add_gatecall!(cct, GateOps.c_u(su, q_tgt, q_ctrl[end]) )
        append!(cct, apply_ncu(q_ctrl[1:end-1], q_aux, q_ctrl[end], GateOps.GateLabel(:x)) )
        Circuit.add_gatecall!(cct, GateOps.c_u(asu, q_tgt, q_ctrl[end]) )
        append!(cct, apply_ncu(q_ctrl[1:end-1], q_aux, q_ctrl[end], GateOps.GateLabel(:x)) )
        append!(cct, apply_ncu(q_ctrl[1:end-1], q_aux, q_tgt, su) )
=#
        apply_cu!(cct, q_ctrl[end], q_tgt, nothing, su)
        append!(cct, apply_ncu(q_ctrl[1:end-1], q_aux, q_ctrl[end], GateOps.GateLabel(:x)) )
        apply_cu!(cct, q_ctrl[end], q_tgt, nothing, asu)
        append!(cct, apply_ncu(q_ctrl[1:end-1], q_aux, q_ctrl[end], GateOps.GateLabel(:x)) )
        append!(cct, apply_ncu(q_ctrl[1:end-1], q_aux, q_tgt, su) )

    end

    return cct

end

# =========================================================================== #
# 
# =========================================================================== #

"""
    apply_ncu!(circuit::Circuit.Circ, q_ctrl::Vector, q_aux::Vector, q_tgt, U::GateOps.GateLabel)

Apply an n-qubit controlled gate operation on the given target qubit.
Ensure the gate corresponding with symbol U is registered with g_map before use.
Appends the GateCall operations to circuit

# Arguments
- `circuit::Circuit.Circ`
- `ctrls::Vector`: 
- `aux::Vector`: 
- `tgt::Int`:
- `U::GateOps.GateLabel`:
"""
function apply_ncu!(circuit::Circuit.Circ, q_ctrl::Vector, q_aux::Vector, q_tgt, U::GateOps.GateLabel)

    # Check global circuit cache for gate
    if ~haskey(Circuit.gate_cache, U )
        error("Gate $(U) does not exist")
    end

    su,asu = get_intermed_gate(U)
    push!(circuit.gate_set, su)
    push!(circuit.gate_set, asu)

    if length(q_ctrl) == 2
        apply_cu!(circuit, q_ctrl[2], q_tgt, nothing, su)
        apply_cx!(circuit, q_ctrl[1], q_ctrl[2], nothing)
        apply_cu!(circuit, q_ctrl[2], q_tgt,  nothing, asu)
        apply_cx!(circuit, q_ctrl[1], q_ctrl[2], nothing)
        apply_cu!(circuit, q_ctrl[1], q_tgt,  nothing, su)

    elseif length(q_ctrl) == 3
        ssu,assu = get_intermed_gate(su)
        push!(circuit.gate_set, ssu)
        push!(circuit.gate_set, assu)

        apply_cu!(circuit, q_ctrl[1], q_tgt, nothing, ssu)
        apply_cx!(circuit, q_ctrl[1], q_ctrl[2], nothing)
        apply_cu!(circuit, q_ctrl[2], q_tgt, nothing, assu)
        apply_cx!(circuit, q_ctrl[1], q_ctrl[2], nothing)
        apply_cu!(circuit, q_ctrl[2], q_tgt, nothing, ssu)
        apply_cx!(circuit, q_ctrl[2], q_ctrl[3], nothing)
        apply_cu!(circuit, q_ctrl[3], q_tgt, nothing, assu)
        apply_cx!(circuit ,q_ctrl[1], q_ctrl[3], nothing)
        apply_cu!(circuit, q_ctrl[3], q_tgt, nothing, ssu)
        apply_cx!(circuit, q_ctrl[2], q_ctrl[3], nothing)
        apply_cu!(circuit, q_ctrl[3], q_tgt, nothing, assu)
        apply_cx!(circuit, q_ctrl[1], q_ctrl[3], nothing)
        apply_cu!(circuit, q_ctrl[3], q_tgt, nothing, ssu)
        
    elseif (length(q_ctrl)>=5) && (length(q_aux)>=length(q_ctrl)-2) && (U == GateOps.GateLabel(:x))
        apply_ncu!(circuit, [q_ctrl[end], q_aux[ 1 + length(q_ctrl)-3 ]], Int[], q_tgt, U)
        for i in reverse(2:length(q_ctrl)-2)
            apply_ncu!(circuit, [q_ctrl[1 + i], q_aux[1 + (i-2)]], Int[], q_aux[1 + (i-1)], U)
        end
        
        apply_ncu!(circuit, [q_ctrl[1], q_ctrl[2]], Int[], q_aux[1], U)
        for i in 2:length(q_ctrl)-2
            apply_ncu!(circuit, [q_ctrl[1 + i], q_aux[1 + (i-2)]], Int[], q_aux[1 + (i-1)], U)
        end
        
        apply_ncu!(circuit, [q_ctrl[end], q_aux[1 + length(q_ctrl) - 3]], Int[], q_tgt, U)
        for i in reverse(2:length(q_ctrl)-2)
            apply_ncu!(circuit, [q_ctrl[1 + i], q_aux[1 + (i-2)]], Int[], q_aux[1 + (i-1)], U)
        end
        
        apply_ncu!(circuit, [q_ctrl[1], q_ctrl[2]], Int[], q_aux[1], U)
        for i in 2:length(q_ctrl)-2
            apply_ncu!(circuit, [q_ctrl[1 + i], q_aux[1 + (i-2)]], Int[], q_aux[1 + (i-1)], U)
        end
    else
        apply_cu!(circuit, q_ctrl[end], q_tgt, nothing, su)
        apply_ncu!(circuit, q_ctrl[1:end-1], q_aux, q_ctrl[end], GateOps.GateLabel(:x))
        apply_cu!(circuit, q_ctrl[end], q_tgt, nothing, asu)
        apply_ncu!(circuit, q_ctrl[1:end-1], q_aux, q_ctrl[end], GateOps.GateLabel(:x))
        apply_ncu!(circuit, q_ctrl[1:end-1], q_aux, q_tgt, su)

    end

    return circuit

end



end
