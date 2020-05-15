module Diffusion

import .NCU
import .Oracle

function apply_h(q_tgt)
    println("h $(q_tgt);")
end
function apply_x(q_tgt)
    println("x $(q_tgt);")
end

"""
    apply_diffusion(q_reg::String, ctrl_indices::Array{Unsigned,1}, tgt_index)

Application of the Grover diffusion operator to marked register.
"""
function apply_diffusion(ctrl_indices::Vector, tgt_index)
    aux_idx = typeof(ctrl_indices)()

    for ctrl in vcat(ctrl_indices, tgt_index)
        apply_h(ctrl)
        apply_x(ctrl)
    end

    NCU.apply_ncu(ctrl_indices, aux_idx, tgt_index, :z)

    for ctrl in vcat(ctrl_indices, tgt_index)
        apply_x(ctrl)
        apply_h(ctrl)
    end
end

end