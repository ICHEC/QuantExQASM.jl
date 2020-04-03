module Decomposition

import ..GateOps
import ..Utils

using PyCall
using LRUCache
using Memoize

# Import decomposition and gate construction utilities from Qiskit, initialised with module __init__ function
const euler_decomposer = PyNULL() 
const u3 = PyNULL()
function __init__()
    copy!(euler_decomposer, pyimport("qiskit.quantum_info.synthesis.one_qubit_decompose").OneQubitEulerDecomposer("U3"))
    copy!(u3, pyimport("qiskit.extensions.standard.u3").U3Gate)
end

"""
    gate_root_adj(angles::GateOps.Gate, decomp_depth::Int)

Calculate intermediate rooted gates using U3 gate definition
@memoize Dict function gate_root_adj( gate::GateOps.Gate, decomp_depth::Int )
    result = sqrt_gate(gate)
    if decomp_depth > 1
        return gate_root_adj(result, decomp_depth-1)
    elseif decomp_depth == 1
        return (result, adjoint_gate(result))
    else
        return (gate, adjoint_gate(gate))
    end
end
"""

"""
    gate_root_adj(angles::GateOps.Gate, decomp_depth::Int)

Calculate intermediate rooted gates using U3 gate definition
"""
@memoize Dict function gate_root_adj( gate::GateOps.Gate )
    result = sqrt_gate(gate)
    return (result, adjoint_gate(result))
end



"""
    matrix_to_u3(mat::Matrix{<:Number} )

Convert SU2 matrix to Euler angles 
"""
function matrix_to_u3(mat::Matrix{<:Number})
    θ, ϕ, λ, phase = euler_decomposer.angles_and_phase(mat)
    return θ, ϕ, λ, phase
end


"""
    u3_to_gate(angles::GateOps.bloch_angles)

Convert Euler angles to SU2 matrix
"""
function u3_to_gate(angles::GateOps.bloch_angles, label::String)
    m = u3(angles.θ, angles.ϕ, angles.λ).to_matrix().*exp(1im*angles.g_phase)
    g = GateOps.Gate(angles, m, label)
    return g
end

"""
    sqrt_gate(angles::GateOps.bloch_angles)

Calculate sqrt of Gate and return as GateOps.Gate
"""
function sqrt_gate(gate::GateOps.Gate)
    mat = u3(gate.angles.θ, gate.angles.ϕ, gate.angles.λ).to_matrix().*exp(1im*gate.angles.g_phase)
    sq_mat = sqrt(mat)
    θ, ϕ, λ, phase = matrix_to_u3(sq_mat)
    return GateOps.Gate(GateOps.bloch_angles(θ, ϕ, λ, phase), sq_mat, "sqrt(" * gate.label * ")")
end

"""
    adjoint_gate(angles::GateOps.bloch_angles)

Calculate adjoint of Gate and return as GateOps.Gate
"""
function adjoint_gate(gate::GateOps.Gate)
    mat = collect(adjoint(u3(gate.angles.θ, gate.angles.ϕ, gate.angles.λ).to_matrix().*exp(1im*gate.angles.g_phase)))
    θ, ϕ, λ, phase = matrix_to_u3(mat)
    return GateOps.Gate(GateOps.bloch_angles(θ, ϕ, λ, phase), mat, "adj(" * gate.label * ")")
end

end