module NCU
using .GateOps

Base.:sqrt(U::Symbol) = Symbol("s"*string(U))
Base.:adjoint(U::Symbol) = Symbol("a"*string(U))

# Start with default gate-set
g_map = Dict{Symbol, Matrix{<:Number}}(
    :x=>([0 1; 1 0].+0im),
    :y=>([0 -1im; 1im 0]),
    :z=>([1 0; 0 -1].+0im),
    :h=>((1/sqrt(2)).*[1 1; 1 -1].+0im)
)

function init_intermed_gates(num_ctrl::Union{Nothing, Int})
    for k in [:z,:y, :x, :h]
        gen_intermed_gates(num_ctrl == nothing ? 8 : num_ctrl, k)
    end
end

function register_gate(U::Symbol, gate::Matrix{<:Number})
    # Some characters are reserved for internal gates
    @assert !(String(U) in ["x","y","z","h","s","t","a","c"])
    g_map[U] = gate
end

function gen_intermed_gates(ctrl_depth::Int, U::Symbol)
    su,asu = get_intermed_gate(U)
    for i in 2:ctrl_depth-1
        su,asu = get_intermed_gate(su)
    end
end

function get_intermed_gate(U::Symbol)
    su = Symbol("s" * String(U))
    asu = Symbol("as" * String(U))

    if haskey(g_map, su)
        SU = g_map[su]
        ASU = g_map[asu]
    else
        SU = sqrt(g_map[U])
        ASU = collect(adjoint(sqrt(g_map[U])))

        # Cache the gates
        g_map[su] = SU
        g_map[asu] = ASU

        # Print the newly generated gates to qasm file for use.
        # Gates can be defined anywhere in file.
        create_gates_nonparam(su, 1)
        create_gates_nonparam(asu, 1)
        create_gates_nonparam(su, 2)
        create_gates_nonparam(asu, 2)
    end
    return su,asu
end

"""
    apply_ncu(ctrls::Vector{Int}, aux::Vector{Int}, tgt::Int, U::Symbol)

Apply an n-qubit controlled gate operation on the given target qubit.
Ensure the gate corresponding with symbol U is registered with g_map before use.

# Arguments
- `ctrls::Vector{Int}`: 
- `aux::Vector{Int}`: 
- `tgt::Int`:
- `U::Symbol`:
"""
function apply_ncu(q_ctrl::Vector, q_aux::Vector, q_tgt, U::Symbol)

    if ~haskey(g_map, U)
        error("Gate $(U) does not exist")
    end

    su,asu = get_intermed_gate(U)

    if length(q_ctrl) == 2
        apply_cu(q_ctrl[2], q_tgt, su)
        apply_cx(q_ctrl[1], q_ctrl[2])
        apply_cu(q_ctrl[2], q_tgt, asu)
        apply_cx(q_ctrl[1], q_ctrl[2])
        apply_cu(q_ctrl[1], q_tgt, su)

    elseif length(q_ctrl) == 3
        #ssu = sqrt(su)
        #assu = adjoint(ssu)
        ssu,assu = get_intermed_gate(su)

        apply_cu(q_ctrl[1], q_tgt, ssu)
        apply_cx(q_ctrl[1], q_ctrl[2])
        apply_cu(q_ctrl[2], q_tgt, assu)
        apply_cx(q_ctrl[1], q_ctrl[2])
        apply_cu(q_ctrl[2], q_tgt, ssu)
        apply_cx(q_ctrl[2], q_ctrl[3])
        apply_cu(q_ctrl[3], q_tgt, assu)
        apply_cx(q_ctrl[1], q_ctrl[3])
        apply_cu(q_ctrl[3], q_tgt, ssu)
        apply_cx(q_ctrl[2], q_ctrl[3])
        apply_cu(q_ctrl[3], q_tgt, assu)
        apply_cx(q_ctrl[1], q_ctrl[3])
        apply_cu(q_ctrl[3], q_tgt, ssu)

    elseif (length(q_ctrl)>=5) && (length(q_aux)>=length(q_ctrl)-2) && (U == :x)
        apply_ncu([q_ctrl[end], q_aux[ 1 + length(q_ctrl)-3 ]], Int[], q_tgt, U) 
        for i in reverse(2:length(q_ctrl)-2)
            apply_ncu([q_ctrl[1 + i], q_aux[1 + (i-2)]], Int[], q_aux[1 + (i-1)], U)
        end
        
        apply_ncu([q_ctrl[1], q_ctrl[2]], Int[], q_aux[1], U)
        for i in 2:length(q_ctrl)-2
            apply_ncu([q_ctrl[1 + i], q_aux[1 + (i-2)]], Int[], q_aux[1 + (i-1)], U)
        end
        
        apply_ncu([q_ctrl[end], q_aux[1 + length(q_ctrl) - 3]], Int[], q_tgt, U)
        for i in reverse(2:length(q_ctrl)-2)
            apply_ncu([q_ctrl[1 + i], q_aux[1 + (i-2)]], Int[], q_aux[1 + (i-1)], U)
        end
        
        apply_ncu([q_ctrl[1], q_ctrl[2]], Int[], q_aux[1], U)
        for i in 2:length(q_ctrl)-2
            apply_ncu([q_ctrl[1 + i], q_aux[1 + (i-2)]], Int[], q_aux[1 + (i-1)], U)
        end

    else
        apply_cu(q_ctrl[end], q_tgt, su)
        apply_ncu(q_ctrl[1:end-1], q_aux, q_ctrl[end], :x)
        apply_cu(q_ctrl[end], q_tgt, asu)
        apply_ncu(q_ctrl[1:end-1], q_aux, q_ctrl[end], :x)
        apply_ncu(q_ctrl[1:end-1], q_aux, q_tgt, su)
    end

end

end
