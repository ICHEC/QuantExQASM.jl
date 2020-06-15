module Oracle

using QuantExQASM.GateOps
using QuantExQASM.NCU
using QuantExQASM.Circuit

function apply_x!(c::Circuit.Circ, tgt, reg::Union{String, Nothing}=nothing) 
    Circuit.add_gatecall!(c, GateOps.pauli_x(tgt, reg))
end

"""
    bitstring_ncu(bitstring::Unsigned, ctrl_indices::Vector, tgt_idx, U::Symbol)

Takes bitstring as the binary pattern and indices as the qubits to operate upon. Applies the appropriate PauliX gates to the control lines to call the NCU with the given matrix 
"""
function bitstring_ncu(cct::Circuit.Circ, bitstring::Integer, ctrl_indices::Vector, tgt_idx, U::GateOps.GateLabel)
    bitmask =  0x1
    aux_idx = typeof(ctrl_indices)()

    # Filter qubit values to mark specified pattern
    for idx in collect(0:length(ctrl_indices))
        if ~( (bitstring & (bitmask << idx) ) != 0)
            if idx < length(ctrl_indices)
                apply_x!(cct, ctrl_indices[idx+1])
            else
                apply_x!(cct, tgt_idx)
            end
        end
    end
    NCU.apply_ncu!(cct, ctrl_indices, aux_idx, tgt_idx, U)
    for idx in collect(0:length(ctrl_indices))
        if ~( (bitstring & (bitmask << idx) ) != 0)
            if idx < length(ctrl_indices)
                apply_x!(cct, ctrl_indices[idx+1])
            else
                apply_x!(cct, tgt_idx)
            end
        end
    end
end

"""
    bitstring_phase_oracle(bitstring::Unsigned, ctrl_indices::Vector, tgt_idx::Unsigned)

Applies PauliX gates to the appropriate lines in the circuit, then applies a n-controlled PauliZ to mark the state.
"""
function bitstring_phase_oracle(cct::Circuit.Circ, bitstring::Integer, ctrl_indices::Vector, tgt_idx)
    return bitstring_ncu(cct, bitstring, ctrl_indices, tgt_idx, GateOps.GateLabel(:z) )
end

end