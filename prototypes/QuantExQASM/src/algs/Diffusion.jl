module Diffusion

using QuantExQASM.NCU
using QuantExQASM.Oracle
using QuantExQASM.Circuit
using QuantExQASM.GateOps

function apply_h!(cct::Circuit.Circ, tgt, reg::Union{String, Nothing}=nothing)
    Circuit.add_gatecall!(cct, GateOps.hadamard(tgt, reg))
end


function apply_x!(c::Circuit.Circ, tgt, reg::Union{String, Nothing}=nothing) 
    Circuit.add_gatecall!(c, GateOps.pauli_x(tgt, reg))
end

"""
    apply_diffusion(q_reg::String, ctrl_indices::Array{Unsigned,1}, tgt_index)

Application of the Grover diffusion operator to marked register.
"""
function apply_diffusion(cct::Circuit.Circ, ctrl_indices::Vector, tgt_index, reg::Union{String, Nothing}=nothing)
    aux_idx = typeof(ctrl_indices)()

    for ctrl in vcat(ctrl_indices, tgt_index)
        apply_h!(cct, ctrl, reg)
        apply_x!(cct, ctrl, reg)
    end

    NCU.apply_ncu!(cct, ctrl_indices, aux_idx, tgt_index, GateOps.GateLabel(:z))

    for ctrl in vcat(ctrl_indices, tgt_index)
        apply_x!(cct, ctrl, reg)
        apply_h!(cct, ctrl, reg)
    end
end

end