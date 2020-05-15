module NCU

function apply_cx(q_ctrl::Int, q_tgt::Int)
    println("X $(q_ctrl) $(q_tgt)")
end

function apply_cu(q_ctrl::Int, q_tgt::Int, U)
    println("$(U) $(q_ctrl) $(q_tgt)")
end

Base.:sqrt(U::Symbol) = Symbol("s"*string(U))
Base.:adjoint(U::Symbol) = Symbol("a"*string(U))
#X = [0 1; 1 0]
X = :x

function apply_ncu(q_ctrl::Vector{Int}, q_aux::Vector{Int}, q_tgt::Int, U)
    
    su = sqrt(U)
    asu = adjoint(su)

    if length(q_ctrl) == 2
        apply_cu(q_ctrl[2], q_tgt, su)
        apply_cx(q_ctrl[1], q_ctrl[2])
        apply_cu(q_ctrl[2], q_tgt, asu)
        apply_cx(q_ctrl[1], q_ctrl[2])
        apply_cu(q_ctrl[1], q_tgt, su)

    elseif length(q_ctrl) == 3
        ssu = sqrt(su)
        assu = adjoint(ssu)

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
        apply_ncu(q_ctrl[1:end-1], q_aux, q_ctrl[end], X)
        apply_cu(q_ctrl[end], q_tgt, asu)
        apply_ncu(q_ctrl[1:end-1], q_aux, q_ctrl[end], X)
        apply_ncu(q_ctrl[1:end-1], q_aux, q_tgt, su)

    end

end

end
