module Grover

import .Oracle
import .Diffusion
import .GateOps

export run_grover

function apply_h(q_tgt)
    println("h $(q_tgt);")
end

function calc_iterations(num_states::Integer)
    return pi*sqrt(num_states)/4
end

function state_init(qubit_indices::Vector)
    for i in qubit_indices
        apply_h(i)
    end
end

function apply_grover_iteration(qubit_indices::Vector)
    return Diffusion.apply_diffusion(qubit_indices[1:end-1], qubit_indices[end])
end

function mark_state(state::Integer, qubit_indices::Vector)
    return Oracle.bitstring_phase_oracle(state, qubit_indices[1:end-1], qubit_indices[end])
end

function run_grover(qubit_indices::Vector, state::Integer)
    num_states = 2^(length(qubit_indices))
    state_init(qubit_indices)
    for i in range(1,stop=ceil(Integer, calc_iterations(num_states)))
        mark_state(state, qubit_indices)
        apply_grover_iteration(qubit_indices)
    end
end

end