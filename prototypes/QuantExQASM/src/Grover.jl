module Grover

import ..Oracle
import ..Diffusion
import ..GateOps

export run_grover

function calc_iterations(num_states::Integer)
    return pi*sqrt(num_states)/4
end

function state_init(q_reg::String, qubit_indices::Vector{<:Integer})
    cct = ""
    for i in qubit_indices
        @debug "H( $(i) ) func=state_init"
        cct *= GateOps.apply_gate_h(q_reg, i)
    end
    return cct
end

function apply_grover_iteration(q_reg::String, qubit_indices::Vector{<:Integer})
    return Diffusion.apply_diffusion(q_reg, qubit_indices[1:end-1], qubit_indices[end])
end

function mark_state(q_reg::String, state::Integer, qubit_indices::Vector{<:Integer})
    return Oracle.bitstring_phase_oracle(q_reg, state, qubit_indices[1:end-1], qubit_indices[end])
end

function run_grover(q_reg::String, qubit_indices::Vector{<:Integer}, state::Integer)
    num_states = 2^(length(qubit_indices))
    result = state_init(q_reg, qubit_indices)
    for i in range(1,stop=ceil(Integer, calc_iterations(num_states)))
        result *= mark_state(q_reg, state, qubit_indices)
        result *= apply_grover_iteration(q_reg, qubit_indices)
    end
    return result
end

end