module Diffusion

import ..GateOps
import ..NCU
import ..Oracle
using DataStructures

"""
    apply_diffusion(q_reg::String, ctrl_indices::Array{Unsigned,1}, tgt_index)

Application of the Grover diffusion operator to marked register.
"""
function apply_diffusion(q_reg::String, ctrl_indices::Vector{<:Integer}, tgt_index)
    result = ""
    uncompute = Stack{String}()
    aux_idx = typeof(ctrl_indices)()

    for ctrl in vcat(ctrl_indices, tgt_index)
        @debug "H( $(ctrl) ) func=apply_diffusion"
        cct = GateOps.apply_gate_h(q_reg, ctrl)

        result *= cct
        push!(uncompute, cct)

        @debug "X( $(ctrl) ) func=apply_diffusion"
        cct = GateOps.apply_gate_x(q_reg, ctrl)
        result *= cct
        push!(uncompute, cct)
    end

    result *= NCU.apply_ncu(q_reg, ctrl_indices, aux_idx, tgt_index, GateOps.default_gates["Z"], 0)

    for ctrl in vcat(ctrl_indices, tgt_index)
        @debug "X( $(ctrl) ) func=apply_diffusion"
        cct = GateOps.apply_gate_x(q_reg, ctrl)
        result *= cct

        @debug "H( $(ctrl) ) func=apply_diffusion"
        cct = GateOps.apply_gate_h(q_reg, ctrl)
        result *= cct
    end

    return result
end

end