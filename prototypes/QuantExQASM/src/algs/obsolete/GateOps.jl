module GateOps

export GateOps
export pauli_x, pauli_y, pauli_z, hadamard, u
export r_x, r_y, r_z, r_phase, swap
export c_pauli_x, c_pauli_y, c_pauli_z, c_u
export c_r_x, c_r_y, c_r_z, c_r_phase

using LinearAlgebra
using Memoize

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

abstract type AGate end
abstract type AGateCall <: AGate end
abstract type AGateLabel <: AGate end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#                       Lightweight gate calls
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

struct GateLabelP{IType<:Integer, NType<:Number} <: AGateLabel
    gate_label::Symbol
    params::Union{Nothing, Dict{String, Union{IType, NType, Bool}}}
    GateLabelP{IType,NType}(gate_label) where {IType<:Integer, NType<:Number} = new(gate_label, nothing)
    GateLabelP{IType,NType}(gate_label, params) where {IType<:Integer, NType<:Number} = new(gate_label, params)
end

struct GateCall1P{IType<:Integer, NType<:Number} <: AGateCall
    gate_label::GateLabelP{IType, NType}
    target::IType
    GateCall1P{IType,NType}(gate_label, target) where {IType<:Integer, NType<:Number} = new(gate_label, target)
end

struct GateCall2P{IType<:Integer, NType<:Number} <: AGateCall
    gate_label::GateLabelP{IType, NType}
    ctrl::IType
    target::IType
    GateCall2P{IType,NType}(gate_label, ctrl, target) where {IType<:Integer, NType<:Number} = new(gate_label, ctrl, target)
end

struct GateCallNP{IType<:Integer, NType<:Number} <: AGateCall
    gate_label::GateLabelP{IType, NType}
    ctrl::Vector{IType}
    target::IType
    GateCallNP{IType,NType}(gate_label, ctrl, target) where {IType<:Integer, NType<:Number} = new(gate_label, ctrl, target)
end

# Useful aliases
GateLabel = GateLabelP{Int64, Float64}
GateCall1 = GateCall1P{Int64, Float64}
GateCall2 = GateCall2P{Int64, Float64}
GateCallN = GateCallNP{Int64, Float64}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#                           Gate implementations
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

struct Gate <: AGateCall
    mat::Matrix{<:Number}
    Gate(mat) = new(mat)
end

struct GateEulerP{NType<:Number} <: AGateCall
    θ::NType
    ϕ::NType
    λ::NType
    phase::NType
    GateEulerP{NType}(θ, ϕ, λ) where {NType<:Number} = new(θ, ϕ, λ, 0.0)
    GateEulerP{NType}(θ, ϕ, λ, phase) where {NType<:Number} = new(θ, ϕ, λ, phase)
    function GateEulerP{NType}(params::Dict{String,NType}) where {NType<:Number}
        θ, ϕ, λ, phase = params["θ"], params["ϕ"], params["λ"], params["phase"]
        return GateEulerP{NType}(θ, ϕ, λ, phase)
    end
end
GateEuler = GateEulerP{Float64}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#                            1Q Gate calls ops
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

function pauli_x(q_target::Int)
    return GateCall1(GateLabel(:x), q_target)
end
function pauli_y(q_target::Int)
    return GateCall1(GateLabel(:y), q_target)
end
function pauli_z(q_target::Int)
    return GateCall1(GateLabel(:z), q_target)
end
function hadamard(q_target::Int)
    return GateCall1(GateLabel(:h), q_target)
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

struct u_gate_schema{NType<:Number}
    schema::Union{Tuple{NType,NType,NType,NType}, Tuple{NType, Bool}}
    type::String
    u_gate_schema{NType}(schema) where {NType<:Number} = length(schema) == 2 ? new(schema, "r") : new(schema, "e")
end

# Arbitrary U matrix can be defined with generating parameters
# Gate schema can be used to indicate construction params 
# ("r" uses given gate, sqrt depth and adjoint; "e" uses euler angles and phase)
function u(label::String, q_target::Int, params::Dict)
    return GateCall1(GateLabel(Symbol(label), params), q_target)
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

function r_x(q_target::Int, theta::Real)
    return GateCall1(GateLabel(Symbol("r_x" * "_" * "angle=" * string(theta)), Dict("angle"=>theta)), q_target )
end
function r_y(q_target::Int, theta::Real)
    return GateCall1(GateLabel(Symbol("r_y" * "_" * "angle=" * string(theta)), Dict("angle"=>theta)), q_target )
end
function r_z(q_target::Int, theta::Real)
    return GateCall1(GateLabel(Symbol("r_z" * "_" * "angle=" * string(theta)), Dict("angle"=>theta)), q_target )
end
function r_phase(q_target::Int, theta::Real)
    return GateCall1(GateLabel(Symbol("r_phase" * "_" * "angle=" * string(theta)), Dict("angle"=>theta)), q_target )
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#                            2Q Gate calls ops
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

function swap(q_target::Int, q_ctrl::Int)
    return GateCall2( GateLabel(:swap), q_target, q_ctrl)
end
function c_pauli_x(q_target::Int, q_ctrl::Int)
    return GateCall2( GateLabel(:c_x), q_target, q_ctrl)
end
function c_pauli_y(q_target::Int, q_ctrl::Int)
    return GateCall2( GateLabel(:c_y), q_target, q_ctrl)
end
function c_pauli_z(q_target::Int, q_ctrl::Int)
    return GateCall2( GateLabel(:c_z), q_target, q_ctrl)
end
function c_u(label::String, q_target::Int, q_ctrl::Int, params::Dict)
    return GateCall2( GateLabel( Symbol(label), params), q_target, q_ctrl)
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

function c_r_x(q_target::Int, q_ctrl::Int, theta::Real)
    return GateCall2( GateLabel( Symbol("c_r_x" * "_" * "angle=" * string(theta)), Dict("angle"=>theta) ), q_target, q_ctrl )
end
function c_r_y(q_target::Int, q_ctrl::Int, theta::Real)
    return GateCall2( GateLabel( Symbol("c_r_y" * "_" * "angle=" * string(theta)), Dict("angle"=>theta) ), q_target, q_ctrl, Dict("theta"=>theta) )
end
function c_r_z(q_target::Int, q_ctrl::Int, theta::Real)
    return GateCall2( GateLabel( Symbol("c_r_z" * "_" * "angle=" * string(theta)), Dict("angle"=>theta) ), q_target, q_ctrl, Dict("theta"=>theta) )
end
function c_r_phase(q_target::Int, q_ctrl::Int, theta::Real)
    return GateCall2( GateLabel( Symbol("c_r_phase" * "_" * "angle=" * string(theta)), Dict("angle"=>theta) ), q_target, q_ctrl, Dict("theta"=>theta) )
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#                           Utility functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

function create_label_r(params::Dict)
    label = params["label"]
    sqrt_depth = params["depth"]
    adj = params["adj"]
    return "$(label)_$(sqrt_depth)_$(adj)"
end
function create_label_e(params::Dict)
    label = params["label"]
    θ = params["θ"]
    ϕ = params["ϕ"]
    λ = params["λ"]
    phase = params["phase"]
    return "$(label)_$(θ)_$(λ)_$(phase)"
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

@memoize function Base.:sqrt(g::Gate)
    return Gate(sqrt(g.mat))
end

@memoize function Base.:adjoint(gate::Gate)
    return Gate(adjoint(gate.mat))
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

"""
Adapated from qiskit's one_qubit_decompose ZYZ method
"""
@memoize function mat_to_euler(gate::Gate)
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
"""
Adapated from qiskit's one_qubit_decompose ZYZ method
"""
@memoize function mat_to_euler(gate::Matrix{<:Number})
    d = 1/sqrt( det( gate ))
    phase = -angle(d)

    su_mat = d*gate
    
    θ = 2.0 * atan( abs(su_mat[2, 1]), abs(su_mat[1, 1]) )
    ϕpλ = 2.0 * ( angle(su_mat[2, 2]) )
    ϕmλ = 2.0 * ( angle(su_mat[2, 1]) )
    ϕ = ( ϕpλ + ϕmλ ) / 2.0
    λ = ( ϕpλ - ϕmλ ) / 2.0

    return (θ ,ϕ, λ, phase)
end


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

```math
\renewcommand{\arraystretch}{1.25}
\begin{pmatrix}
\cos\frac{\theta}{2} & e^{-i\lambda}\sin\frac{\theta}{2} \\
e^{i\phi}\sin\frac{\theta}{2} & e^{i(\phi + \lambda)}\cos\frac{\theta}{2} \\
\end{pmatrix}
```
@memoize function euler_to_mat(θ::Number, ϕ::Number, λ::Number, phase::Number)
    a =  exp( 1im*phase ) * cos(θ/2)
    b =  exp(-1im*(λ-phase) ) * sin(θ/2)
    c =  exp( 1im*(ϕ+phase) ) * sin(θ/2)
    d =  exp( 1im*( ϕ + λ + phase) ) * cos(θ/2)
    return [a b; c d]
end
@memoize function euler_to_mat(params::GateEulerP)
    a =  exp( 1im*params.phase ) * cos(params.θ/2)
    b =  exp(-1im*( params.λ-params.phase) ) * sin(params.θ/2)
    c =  exp( 1im*( params.ϕ+params.phase) ) * sin(params.θ/2)
    d =  exp( 1im*( params.ϕ + params.λ + params.phase) ) * cos(params.θ/2)
    return [a b; c d]
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

end



@memoize function mat_to_euler(unitary)
    """Finds rotation angles (a,b,c,d) in the decomposition u=exp(id)*Rz(c).Ry(b).Rz(a).
    Ported to Julia from qiskit.
    """
    u00 = unitary[1, 1]
    u01 = unitary[1, 2]
    u10 = unitary[2, 1]
    u11 = unitary[2, 2]

    tol = 1e-8

    # Handle special case if the entry (0,0) of the unitary is equal to zero
    if abs(u00) <  tol
        # Note that u10 can't be zero, since u is unitary (and u00 == 0)
        gamma = angle(-u01 / u10)
        delta = angle(u01 * exp(-1im * gamma / 2))
        return 0., pi, -gamma, delta
    end
    # Handle special case if the entry (0,1) of the unitary is equal to zero
    if abs(u01) <tol
        # Note that u11 can't be zero, since u is unitary (and u01 == 0)
        gamma = angle(u00 / u11)
        delta = angle(u00 * exp(-1im * gamma / 2))
        return 0., 0., -gamma, delta
    end
    beta = 2 * acos(abs(u00))
    if sin(beta / 2) - cos(beta / 2) > 0
        gamma = angle(-u00 / u10)
        alpha = angle(u00 / u01)
    else
        gamma = -angle(-u10 / u00)
        alpha = -angle(u01 / u00)
    end
    delta = angle(u00 * exp(-1im * (alpha + gamma) / 2))
    # The decomposition works with another convention for the rotation gates
    # (the one using negative angles).
    # Therefore, we have to take the inverse of the angles at the end.
    return -alpha, -beta, -gamma, delta

end

@memoize function euler_to_mat(a,b,c,d)
    """Creates unitary from angles defined by mat_to_euler function as u=exp(id)*Rz(c).Ry(b).Rz(a).
    """
    u00 = exp(-1im*(c+a)/2)*cos(b/2)
    u01 = -exp(1im*(a-c)/2)*sin(b/2)
    u10 = exp(1im*(c-a)/2)*sin(b/2)
    u11 = exp(1im*(a+c)/2)*cos(b/2)
    return [u00 u01; u10 u11].*exp(1im*d)
end

function rz(theta)
    return [exp(-1im*theta/2) 0; 0 exp(1im*theta/2)]
end
function ry(theta)
    return [[cos(theta/2) -sin(theta/2)]; [sin(theta/2) cos(theta/2)]]
end