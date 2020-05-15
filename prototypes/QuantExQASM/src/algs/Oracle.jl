module Oracle

import .NCU

function apply_x(q_tgt)
    println("x $(q_tgt);")
end

"""
    bitstring_ncu(bitstring::Unsigned, ctrl_indices::Vector, tgt_idx, U::Symbol)

Takes bitstring as the binary pattern and indices as the qubits to operate upon. Applies the appropriate PauliX gates to the control lines to call the NCU with the given matrix 
"""
function bitstring_ncu(bitstring::Integer, ctrl_indices::Vector, tgt_idx, U::Symbol)
    bitmask =  0x1
    aux_idx = typeof(ctrl_indices)()

    # Filter qubit values to mark specified pattern
    for idx in collect(0:length(ctrl_indices))
        if ~( (bitstring & (bitmask << idx) ) != 0)
            if idx < length(ctrl_indices)
                apply_x(ctrl_indices[idx+1])
            else
                apply_x(tgt_idx)
            end
        end
    end
    NCU.apply_ncu(ctrl_indices, aux_idx, tgt_idx, U)
    for idx in collect(0:length(ctrl_indices))
        if ~( (bitstring & (bitmask << idx) ) != 0)
            if idx < length(ctrl_indices)
                apply_x(ctrl_indices[idx+1])
            else
                apply_x(tgt_idx)
            end
        end
    end
end

"""
    bitstring_phase_oracle(bitstring::Unsigned, ctrl_indices::Vector, tgt_idx::Unsigned)

Applies PauliX gates to the appropriate lines in the circuit, then applies a n-controlled PauliZ to mark the state.
"""
function bitstring_phase_oracle(bitstring::Integer, ctrl_indices::Vector, tgt_idx)
    return bitstring_ncu(bitstring, ctrl_indices, tgt_idx, :z )
end

end