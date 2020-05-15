module Decomposition

using LinearAlgebra

import ..GateOps
import "../io/Utils"

using Memoize

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

Base.:sqrt(a::GateOps.Gate) = sqrt(a.mat)

"""
Define sqrt for Gate and GateLabel structs
"""
@memoize function Base.:sqrt(g::GateOps.Gate)
    return GateOps.Gate(sqrt(g.mat))
end
@memoize function Base.:sqrt(g::GateOps.GateLabelP)
    m = sqrt(GateOps.GateMap[g])
    GateOps.GateLabelP
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

"""
Adapated from qiskit's one_qubit_decompose ZYZ method
"""
@memoize function mat_to_euler(gate::GateOps.Gate)
    d = 1/sqrt( det( gate.mat ))
    phase = -angle(d)

    su_mat = d*gate.mat
    
    θ = 2.0 * atan( abs(su_mat[2, 1]), abs(su_mat[1, 1]) )
    ϕpλ = 2.0 * ( angle(su_mat[2, 2]) )
    ϕmλ = 2.0 * ( angle(su_mat[2, 1]) )
    ϕ = ( ϕpλ + ϕmλ ) / 2.0
    λ = ( ϕpλ - ϕmλ ) / 2.0

    return (θ ,ϕ, λ, phase)
end

function mat_to_euler(gate::GateOps.GateLabel)
    return mat_to_euler(GateOps.GateMap[gate])
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
"""
```math
\renewcommand{\arraystretch}{1.25}
\begin{pmatrix}
\cos\frac{\theta}{2} & e^{-i\lambda}\sin\frac{\theta}{2} \\
e^{i\phi}\sin\frac{\theta}{2} & e^{i(\phi + \lambda)}\cos\frac{\theta}{2} \\
\end{pmatrix}
```
"""
function euler_to_mat(θ::Number, ϕ::Number, λ::Number, phase::Number)
    a =  exp( 1im*phase ) * cos(θ/2)
    b =  exp(-1im*(λ-phase) ) * sin(θ/2)
    c =  exp( 1im*(ϕ+phase) ) * sin(θ/2)
    d =  exp( 1im*( ϕ + λ + phase) ) * cos(θ/2)
    return [a b; c d]
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

function sqrt_gate(gate::GateOps.Gate)
    return GateOps.Gate(sqrt(gate.mat))
end
function sqrt_gate(gate::GateOps.GateLabel)
    return sqrt( GateOps.GateMap[gate].mat )
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

function adjoint_gate(gate::GateOps.Gate)
    return GateOps.Gate(adjoint(gate.mat))
end
function adjoint_gate(gate::GateOps.GateLabel)
    return adjoint( GateOps.GateMap[gate].mat )
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

end