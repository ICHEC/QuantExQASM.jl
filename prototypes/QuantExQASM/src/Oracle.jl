module Oracle

import ..GateOps
import ..NCU
using DataStructures

"""
    bitstring_ncu(q_reg::String, bitstring::Unsigned, ctrl_indices::Array{Unsigned,1}, tgt_idx::Unsigned, U::GateOps.Gate)

Takes bitstring as the binary pattern and indices as the qubits to operate upon. Applies the appropriate PauliX gates to the control lines to call the NCU with the given matrix 
"""
function bitstring_ncu(q_reg::String, bitstring::Integer, ctrl_indices::Vector{<:Integer}, tgt_idx::Integer, U::GateOps.Gate)
    bitmask =  0x1
    uncompute = Stack{String}()
    result = ""
    aux_idx = typeof(ctrl_indices)()

    # Filter qubit values to mark specified pattern
    for idx in collect(0:length(ctrl_indices))
        if ~( (bitstring & (bitmask << idx) ) != 0)
            if idx < length(ctrl_indices)
                @debug "X( $(ctrl_indices[idx+1]) ) func=bitstring_ncu"
                cct = GateOps.apply_gate_x(q_reg, ctrl_indices[idx+1])
            else
                @debug "X( $(tgt_idx) ) func=bitstring_ncu"
                cct = GateOps.apply_gate_x(q_reg, tgt_idx)
            end
            result *= cct
            push!(uncompute, cct)
        end
    end
    result *= NCU.apply_ncu(q_reg, ctrl_indices, aux_idx, tgt_idx, U, 0)

    for idx in collect(0:length(ctrl_indices))
        if ~( (bitstring & (bitmask << idx) ) != 0)
            if idx < length(ctrl_indices)
                @debug "X( $(ctrl_indices[idx+1]) ) func=bitstring_ncu"
                cct = GateOps.apply_gate_x(q_reg, ctrl_indices[idx+1])
            else
                @debug "X( $(tgt_idx) ) func=bitstring_ncu"
                cct = GateOps.apply_gate_x(q_reg, tgt_idx)
            end
            result *= cct
        end
    end
    return result
end

"""
    bitstring_phase_oracle(q_reg::String, bitstring::Unsigned, ctrl_indices::Array{Unsigned}, tgt_idx::Unsigned)

Applies PauliX gates to the appropriate lines in the circuit, then applies a n-controlled PauliZ to mark the state.
"""
function bitstring_phase_oracle(q_reg::String, bitstring::Integer, ctrl_indices::Vector{<:Integer}, tgt_idx::Integer)
    return bitstring_ncu(q_reg, bitstring, ctrl_indices, tgt_idx, GateOps.default_gates["Z"] )
end

end