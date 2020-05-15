module Grover

using QuantExQASM.Oracle
using QuantExQASM.Diffusion
using QuantExQASM.GateOps
using QuantExQASM.Circuit

export run_grover

function apply_h!(cct::Circuit.Circ, tgt, reg::Union{String, Nothing}=nothing)
    Circuit.add_gatecall!(cct, GateOps.hadamard(tgt, reg))
end

function apply_x!(cct::Circuit.Circ, tgt, reg::Union{String, Nothing}=nothing) 
    Circuit.add_gatecall!(c, GateOps.pauli_x(tgt, reg))
end


function calc_iterations(num_states::Integer)
    return pi*sqrt(num_states)/4
end

function state_init!(cct::Circuit.Circ, qubit_indices::Vector)
    for i in qubit_indices
        apply_h!(cct, i)
    end
end

function apply_grover_iteration!(cct::Circuit.Circ, qubit_indices::Vector)
    return Diffusion.apply_diffusion(cct, qubit_indices[1:end-1], qubit_indices[end])
end

function mark_state!(cct::Circuit.Circ, state::Integer, qubit_indices::Vector)
    return Oracle.bitstring_phase_oracle(cct::Circuit.Circ, state, qubit_indices[1:end-1], qubit_indices[end])
end

function run_grover!(cct::Circuit.Circ, qubit_indices::Vector, state::Integer)
    num_states = 2^(length(qubit_indices))
    state_init!(cct, qubit_indices)
    for i in range(1,stop=ceil(Integer, calc_iterations(num_states)))
        mark_state!(cct, state, qubit_indices)
        apply_grover_iteration!(cct, qubit_indices)
    end
end

end