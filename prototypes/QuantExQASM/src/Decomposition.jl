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
    gate_root_adj(gate::Matrix{<:Number}, decomp_depth::Int)

Calculate intermediate rooted gates and cache intermediates for faster recall using matrix gate definition
"""
@memoize Dict function gate_root_adj(gate::Matrix{<:Number}, decomp_depth::Int)
    result = sqrt(gate)
    if decomp_depth > 1
        return gate_root_adj(result, decomp_depth-1)
    elseif decomp_depth == 1
        return (result, collect(adjoint(result)))
    else
        return (gate, collect(adjoint(gate)))
    end
end


"""
    gate_root_adj(angles::GateOps.bloch_angles, decomp_depth::Int)

Calculate intermediate rooted gates using U3 gate definition
"""
@memoize Dict function gate_root_adj(angles::GateOps.bloch_angles, decomp_depth::Int)
    g_mat = u3_to_matrix(angles)
    g, g_adj = gate_root_adj(g_mat, decomp_depth)
    return (matrix_to_u3(g), matrix_to_u3(g_adj))
end


"""
    u3_to_matrix(angles::GateOps.bloch_angles)

Convert Euler angles to SU2 matrix
"""
function u3_to_matrix(angles::GateOps.bloch_angles)
    return u3(angles.θ, angles.ϕ, angles.λ).to_matrix()
end

"""
    matrix_to_u3(mat::Matrix{<:Number} )

Convert SU2 matrix to Euler angles 
"""
function matrix_to_u3(mat::Matrix{<:Number})
    θ, ϕ, λ = euler_decomposer._angles(mat)
    return GateOps.bloch_angles(θ, ϕ, λ)
end

end