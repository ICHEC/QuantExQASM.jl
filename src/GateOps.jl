module GateOps

export GateOps
export pauli_x, pauli_y, pauli_z, hadamard, u
export r_x, r_y, r_z, r_phase, swap
export c_pauli_x, c_pauli_y, c_pauli_z, c_u
export c_r_x, c_r_y, c_r_z, c_r_phase

using LinearAlgebra
using Memoize

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#                       Redirect I/O streams
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#(rd, wr) = redirect_stdout();

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#                       Lightweight gate calls
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

"Abstract Gate type"
abstract type AGate end
"Abstract Gate Label, for unique gates"
abstract type AGateLabel <: AGate end
"Abstract Gate Call, for tracking Gate labels applied to specific qubits"
abstract type AGateCall <: AGate end

"Parametric Gate label. Tracks the gate symbol (:x,:y,:z, etc), and arbitrary parameters"
struct GateLabelP{IType<:Integer, NType<:Number} <: AGateLabel
    label::Symbol
    params::Union{Nothing, Dict{String, Union{IType, NType, Bool}}}
    GateLabelP{IType,NType}(gate_label) where {IType<:Integer, NType<:Number} = new(gate_label, nothing)
    GateLabelP{IType,NType}(gate_label, params) where {IType<:Integer, NType<:Number} = new(gate_label, params)
end

"Parametric single qubit Gate call. Has a GateLabelP, target qubit index and register label"
struct GateCall1P{IType<:Integer, NType<:Number} <: AGateCall
    gate_label::GateLabelP{IType, NType}
    target::IType
    reg::Union{String, Nothing}
    GateCall1P{IType,NType}(gate_label, target, reg) where {IType<:Integer, NType<:Number} = new(gate_label, target, reg)
    GateCall1P{IType,NType}(gate_label, target) where {IType<:Integer, NType<:Number} = new(gate_label, target, nothing)
end

"""
Parametric two qubit Gate call. 
Has a GateLabelP, control and target qubit indices, register label, and base gate (assumes controlled U).
"""
struct GateCall2P{IType<:Integer, NType<:Number} <: AGateCall
    gate_label::GateLabelP{IType, NType}
    ctrl::IType
    target::IType
    base_gate::Union{GateCall1P{IType, NType}, GateCall2P{IType, NType}, Nothing}
    #base_gate::Union{GateLabelP{IType, NType}, Nothing}
    reg::Union{String, Nothing}
    GateCall2P{IType,NType}(gate_label, ctrl, target, base_gate, reg) where {IType<:Integer, NType<:Number} = 
    new(gate_label, ctrl, target, base_gate, reg)


    GateCall2P{IType,NType}(gate_label, ctrl, target, base_gate) where {IType<:Integer, NType<:Number} = 
    new(gate_label, ctrl, target, base_gate, nothing)


    GateCall2P{IType,NType}(gate_label, ctrl, target) where {IType<:Integer, NType<:Number} = 
    new(gate_label, ctrl, target, nothing, nothing)
end

"""
Parametric n-qubit Gate call. 
Has a GateLabelP, control vector indices, target qubit index, register label, and base gate.
"""
struct GateCallNP{IType<:Integer, NType<:Number} <: AGateCall
    gate_label::GateLabelP{IType, NType}
    ctrl::Vector{IType}
    target::IType
    base_gate::Union{GateCall1P{IType, NType}, GateCall2P{IType, NType}, GateCallNP{IType, NType}, Nothing}
    reg::Union{String,Nothing}
    GateCallNP{IType,NType}(gate_label, ctrl, target, base_gate, reg) where {IType<:Integer, NType<:Number} = new(gate_label, ctrl, target, base_gate{IType,NType}, reg)
    GateCallNP{IType,NType}(gate_label, ctrl, target, base_gate) where {IType<:Integer, NType<:Number} = new(gate_label, ctrl, target, base_gate{IType,NType}, "")
end

# Useful aliases
GateLabel = GateLabelP{Int64, Float64}
GateCall1 = GateCall1P{Int64, Float64}
GateCall2 = GateCall2P{Int64, Float64}
GateCallN = GateCallNP{Int64, Float64}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#                           Gate implementations
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

"Wrap of Matrix to fit type system hierarchy"
struct Gate <: AGateCall
    mat::Matrix{<:Number}
    Gate(mat) = new(mat)
end

"Placeholder struct for Euler-angle gate"
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

"""
    pauli_x(q_target::Int, register::Union{String, Nothing}=nothing)

Generate a single qubit Pauli-x GateCall (GateCall1P) applied to the target qubit (on given register, if provided)

# Examples
```julia-repl
julia> GateOps.pauli_x(0, "q")
QuantExQASM.GateOps.GateCall1P{Int64,Float64}(QuantExQASM.GateOps.GateLabelP{Int64,Float64}(:x, nothing), 0, "q")
```
"""
function pauli_x(q_target::Int, register::Union{String, Nothing}=nothing)
    if register == nothing
        return GateCall1(GateLabel(:x), q_target)
    else
        return GateCall1(GateLabel(:x), q_target, register)
    end
end

"""
    pauli_y(q_target::Int, register::Union{String, Nothing}=nothing)

Generate a single qubit Pauli-y GateCall (GateCall1P) applied to the target qubit (on given register, if provided)

# Examples
```julia-repl
julia> GateOps.pauli_y(0, "q")
QuantExQASM.GateOps.GateCall1P{Int64,Float64}(QuantExQASM.GateOps.GateLabelP{Int64,Float64}(:y, nothing), 0, "q")
```
"""
function pauli_y(q_target::Int, register::Union{String, Nothing}=nothing)
    if register == nothing
        return GateCall1(GateLabel(:y), q_target)
    else
        return GateCall1(GateLabel(:y), q_target, register)
    end
end

"""
    pauli_z(q_target::Int, register::Union{String, Nothing}=nothing)

Generate a single qubit Pauli-z GateCall (GateCall1P) applied to the target qubit (on given register, if provided)

# Examples
```julia-repl
julia> GateOps.pauli_z(0, "q")
QuantExQASM.GateOps.GateCall1P{Int64,Float64}(QuantExQASM.GateOps.GateLabelP{Int64,Float64}(:z, nothing), 0, "q")
```
"""
function pauli_z(q_target::Int, register::Union{String, Nothing}=nothing)
    if register == nothing
        return GateCall1(GateLabel(:z), q_target)
    else
        return GateCall1(GateLabel(:z), q_target, register)
    end
end

"""
    hadamard(q_target::Int, register::Union{String, Nothing}=nothing)

Generate a single qubit Hadamard GateCall (GateCall1P) applied to the target qubit (on given register, if provided)

# Examples
```julia-repl
julia> GateOps.hadamard(0, "q")
GateOps.GateCall1P{Int64,Float64}(GateOps.GateLabelP{Int64,Float64}(:h, nothing), 0, "q")
```
"""
function hadamard(q_target::Int, register::Union{String, Nothing}=nothing)
    if register == nothing
        return GateCall1(GateLabel(:h), q_target)
    else
        return GateCall1(GateLabel(:h), q_target, register)
    end
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

struct u_gate_schema{NType<:Number}
    schema::Union{Tuple{NType,NType,NType,NType}, Tuple{NType, Bool}}
    type::String
    u_gate_schema{NType}(schema) where {NType<:Number} = length(schema) == 2 ? new(schema, "r") : new(schema, "e")
end

"""
    u(label::GateLabel, q_target::Int, register::Union{String, Nothing}=nothing)

Generate a single qubit arbitrary unitary GateCall (GateCall1P) applied to the target qubit (on given register, if provided)

# Examples
```julia-repl
julia> GateOps.u(GateOps.GateLabel(:mygate), 1, "qr")
GateOps.GateCall1P{Int64,Float64}(GateOps.GateLabelP{Int64,Float64}(:mygate, nothing), 1, "qr")
```
"""
function u(label::GateLabel, q_target::Int, register::Union{String, Nothing}=nothing)
    return GateCall1( label, q_target, register)
end
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

"""
    r_x(q_target::Int, theta::Real, register::Union{String, Nothing}=nothing)

Generate a single qubit R_x(θ) GateCall (GateCall1P) applied to the target qubit (on given register, if provided)

# Examples
```julia-repl
julia> GateOps.r_x(3,pi/2)
GateOps.GateCall1P{Int64,Float64}(GateOps.GateLabelP{Int64,Float64}(Symbol("r_x_angle=1.5707963267948966"), Dict{String,Union{Bool, Float64, Int64}}("angle"=>1.5708)), 3, nothing)
```
"""
function r_x(q_target::Int, theta::Number, register::Union{String, Nothing}=nothing)
    if register == nothing
        return GateCall1(GateLabel(Symbol("r_x" * "_" * "angle=" * string(theta)), Dict("angle"=>theta)), q_target )
    else
        return GateCall1(GateLabel(Symbol("r_x" * "_" * "angle=" * string(theta)), Dict("angle"=>theta)), q_target, register )
    end
end

"""
    r_y(q_target::Int, theta::Real, register::Union{String, Nothing}=nothing)

Generate a single qubit R_y(θ) GateCall (GateCall1P) applied to the target qubit (on given register, if provided)

# Examples
```julia-repl
julia> GateOps.r_y(3,pi/2)
GateOps.GateCall1P{Int64,Float64}(GateOps.GateLabelP{Int64,Float64}(Symbol("r_y_angle=1.5707963267948966"), Dict{String,Union{Bool, Float64, Int64}}("angle"=>1.5708)), 3, nothing)
```
"""
function r_y(q_target::Int, theta::Number, register::Union{String, Nothing}=nothing)
    if register == nothing
        return GateCall1(GateLabel(Symbol("r_y" * "_" * "angle=" * string(theta)), Dict("angle"=>theta)), q_target )
    else
        return GateCall1(GateLabel(Symbol("r_y" * "_" * "angle=" * string(theta)), Dict("angle"=>theta)), q_target, register )
    end
end

"""
    r_z(q_target::Int, theta::Real, register::Union{String, Nothing}=nothing)

Generate a single qubit R_z(θ) GateCall (GateCall1P) applied to the target qubit (on given register, if provided)

# Examples
```julia-repl
julia> GateOps.r_z(3,pi/2)
GateOps.GateCall1P{Int64,Float64}(GateOps.GateLabelP{Int64,Float64}(Symbol("r_z_angle=1.5707963267948966"), Dict{String,Union{Bool, Float64, Int64}}("angle"=>1.5708)), 3, nothing)
```
"""
function r_z(q_target::Int, theta::Number, register::Union{String, Nothing}=nothing)
    if register == nothing
        return GateCall1(GateLabel(Symbol("r_z" * "_" * "angle=" * string(theta)), Dict("angle"=>theta)), q_target )
    else
        return GateCall1(GateLabel(Symbol("r_z" * "_" * "angle=" * string(theta)), Dict("angle"=>theta)), q_target, register )
    end
end

"""
    r_phase(q_target::Int, theta::Real, register::Union{String, Nothing}=nothing)

Generate a single qubit phase shift GateCall ( diag(1, exp(1im*theta)) , GateCall1P) applied to the target qubit (on given register, if provided)

# Examples
```julia-repl
julia> GateOps.r_phase(3,pi/2)
GateOps.GateCall1P{Int64,Float64}(GateOps.GateLabelP{Int64,Float64}(Symbol("r_phase_angle=1.5707963267948966"), Dict{String,Union{Bool, Float64, Int64}}("angle"=>1.5708)), 3, nothing)
```
"""
function r_phase(q_target::Int, theta::Number, register::Union{String, Nothing}=nothing)
    if register == nothing
        return GateCall1(GateLabel(Symbol("r_phase" * "_" * "angle=" * string(theta)), Dict("angle"=>theta)), q_target )
    else
        return GateCall1(GateLabel(Symbol("r_phase" * "_" * "angle=" * string(theta)), Dict("angle"=>theta)), q_target, register )
    end
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#                            2Q Gate calls ops
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

function swap(q_target::Int, q_ctrl::Int, register::Union{String, Nothing}=nothing)
    return GateCall2( GateLabel(:swap), q_target, q_ctrl)
end

"""
    c_pauli_x(q_target::Int, q_ctrl::Int, register::Union{String, Nothing}=nothing)

Generate a controlled Pauli-x (two qubit) GateCall (GateCall2P), controlled on the index `q_ctrl` applied to the target `q_target` (on given register, if provided).

# Examples
```julia-repl
julia> GateOps.c_pauli_x(0, 1, "qreg")
GateOps.GateCall2P{Int64,Float64}(GateOps.GateLabelP{Int64,Float64}(:c_x, nothing), 0, 1, GateOps.GateCall1P{Int64,Float64}(GateOps.GateLabelP{Int64,Float64}(:x, nothing), 0, "qreg"), "qreg")
```
"""
function c_pauli_x(q_target::Int, q_ctrl::Int, register::Union{String, Nothing}=nothing)
    return GateCall2( GateLabel(:c_x), q_target, q_ctrl, pauli_x(q_target,register), register)
end

"""
    c_pauli_y(q_target::Int, q_ctrl::Int, register::Union{String, Nothing}=nothing)

Generate a controlled Pauli-y (two qubit) GateCall (GateCall2P), controlled on the index `q_ctrl` applied to the target `q_target` (on given register, if provided).

# Examples
```julia-repl
julia> GateOps.c_pauli_y(0, 1, "qreg")
GateOps.GateCall2P{Int64,Float64}(GateOps.GateLabelP{Int64,Float64}(:c_y, nothing), 0, 1, GateOps.GateCall1P{Int64,Float64}(GateOps.GateLabelP{Int64,Float64}(:y, nothing), 0, "qreg"), "qreg")
```
"""
function c_pauli_y(q_target::Int, q_ctrl::Int, register::Union{String, Nothing}=nothing)
    return GateCall2( GateLabel(:c_y), q_target, q_ctrl, pauli_y(q_target,register), register)
end

"""
    c_pauli_z(q_target::Int, q_ctrl::Int, register::Union{String, Nothing}=nothing)

Generate a controlled Pauli-x (two qubit) GateCall (GateCall2P), controlled on the index `q_ctrl` applied to the target `q_target` (on given register, if provided).

# Examples
```julia-repl
julia> GateOps.c_pauli_z(0, 1, "qreg")
GateOps.GateCall2P{Int64,Float64}(GateOps.GateLabelP{Int64,Float64}(:c_z, nothing), 0, 1, GateOps.GateCall1P{Int64,Float64}(GateOps.GateLabelP{Int64,Float64}(:z, nothing), 0, "qreg"), "qreg")
```
"""
function c_pauli_z(q_target::Int, q_ctrl::Int, register::Union{String, Nothing}=nothing)
    return GateCall2( GateLabel(:c_z), q_target, q_ctrl, pauli_z(q_target,register), register)
end

function c_u(label::String, q_target::Int, q_ctrl::Int, params::Dict, register::Union{String, Nothing}=nothing)
    return GateCall2( GateLabel( Symbol(label), params), q_target, q_ctrl, register)
end

"""
    c_u(label::GateLabel, q_target::Int, q_ctrl::Int, gc::GateCall1, register::Union{String, Nothing}=nothing)

Generate a controlled unitary (two qubit) GateCall (GateCall2P), controlled on the index `q_ctrl` applied to the target `q_target` (on given register, if provided).

# Examples
```julia-repl
julia> GateOps.c_u(GateOps.GateLabel(:myCU), 0, 1, GateOps.pauli_x(0), "q")
GateOps.GateCall2P{Int64,Float64}(GateOps.GateLabelP{Int64,Float64}(:myCU, nothing), 0, 1, GateOps.GateCall1P{Int64,Float64}(GateOps.GateLabelP{Int64,Float64}(:x, nothing), 0, nothing), "q")
```
"""
function c_u(label::GateLabel, q_target::Int, q_ctrl::Int, gc::GateCall1, register::Union{String, Nothing}=nothing)
    return GateCall2( label, q_target, q_ctrl, gc, register)
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

"""
    c_r_x(q_target::Int, q_ctrl::Int, theta::Real, register::Union{String, Nothing}=nothing)

Generate a controlled rotation about x (exp(iθσ_x/2)) GateCall (GateCall2P), controlled on the index `q_ctrl` applied to the target `q_target` (on given register, if provided).

# Examples
```julia-repl
julia> QuantExQASM.GateOps.c_r_x( 0, 1, pi/2, "q")
QuantExQASM.GateOps.GateCall2P{Int64,Float64}(QuantExQASM.GateOps.GateLabelP{Int64,Float64}(Symbol("c_r_x_angle=1.5707963267948966"), Dict{String,Union{Bool, Float64, Int64}}("angle" => 1.5707963267948966)), 0, 1, nothing, nothing)
```
"""
function c_r_x(q_target::Int, q_ctrl::Int, theta::Real, register::Union{String, Nothing}=nothing)
    return GateCall2( GateLabel( Symbol("c_r_x" * "_" * "angle=" * string(theta)), Dict("angle"=>theta) ), q_target, q_ctrl )
end

"""
    c_r_y(q_target::Int, q_ctrl::Int, theta::Real, register::Union{String, Nothing}=nothing)

Generate a controlled rotation about y (exp(iθσ_y/2)) GateCall (GateCall2P), controlled on the index `q_ctrl` applied to the target `q_target` (on given register, if provided).

# Examples
```julia-repl
julia> QuantExQASM.GateOps.c_r_y( 0, 1, pi/3, "q")
QuantExQASM.GateOps.GateCall2P{Int64,Float64}(QuantExQASM.GateOps.GateLabelP{Int64,Float64}(Symbol("c_r_y_angle=1.0471975511965976"), Dict{String,Union{Bool, Float64, Int64}}("angle" => 1.0471975511965976)), 0, 1, nothing, nothing)
```
"""
function c_r_y(q_target::Int, q_ctrl::Int, theta::Real, register::Union{String, Nothing}=nothing)
    return GateCall2( GateLabel( Symbol("c_r_y" * "_" * "angle=" * string(theta)), Dict("angle"=>theta) ), q_target, q_ctrl )
end

"""
    c_r_z(q_target::Int, q_ctrl::Int, theta::Real, register::Union{String, Nothing}=nothing)

Generate a controlled rotation about z (exp(iθσ_z/2)) GateCall (GateCall2P), controlled on the index `q_ctrl` applied to the target `q_target` (on given register, if provided).

# Examples
```julia-repl
julia> QuantExQASM.GateOps.c_r_z( 0, 1, pi/4, "q")
QuantExQASM.GateOps.GateCall2P{Int64,Float64}(QuantExQASM.GateOps.GateLabelP{Int64,Float64}(Symbol("c_r_z_angle=0.7853981633974483"), Dict{String,Union{Bool, Float64, Int64}}("angle" => 0.7853981633974483)), 0, 1, nothing, nothing)
```
"""
function c_r_z(q_target::Int, q_ctrl::Int, theta::Real, register::Union{String, Nothing}=nothing)
    return GateCall2( GateLabel( Symbol("c_r_z" * "_" * "angle=" * string(theta)), Dict("angle"=>theta) ), q_target, q_ctrl )
end

"""
    c_r_phase(q_target::Int, q_ctrl::Int, theta::Real, register::Union{String, Nothing}=nothing)

Generate a controlled phase rotation about (controlled [1 0; 0 exp(iθ)] GateCall (GateCall2P), controlled on the index `q_ctrl` applied to the target `q_target` (on given register, if provided).

# Examples
```julia-repl
julia> QuantExQASM.GateOps.c_r_phase( 0, 1, pi/7, "q")
QuantExQASM.GateOps.GateCall2P{Int64,Float64}(QuantExQASM.GateOps.GateLabelP{Int64,Float64}(Symbol("c_r_phase_angle=0.4487989505128276"), Dict{String,Union{Bool, Float64, Int64}}("angle" => 0.4487989505128276)), 0, 1, nothing, nothing)
```
"""
function c_r_phase(q_target::Int, q_ctrl::Int, theta::Real, register::Union{String, Nothing}=nothing)
    return GateCall2( GateLabel( Symbol("c_r_phase" * "_" * "angle=" * string(theta)), Dict("angle"=>theta) ), q_target, q_ctrl )
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

"""Finds rotation angles (a,b,c,d) in the decomposition u=exp(id)*Rz(c).Ry(b).Rz(a).
Direct port to Julia from qiskit: https://qiskit.org/documentation/_modules/qiskit/extensions/quantum_initializer/squ.html
"""
@memoize function mat_to_euler(unitary::Matrix{<:Number})

    u00 = unitary[1, 1]
    u01 = unitary[1, 2]
    u10 = unitary[2, 1]
    u11 = unitary[2, 2]

    tol = 1e-8 # Tolerance defines proximity to 0; assume 0 if less than tol

    # Handle special case if the entry (0,0) of the unitary is equal to zero
    if abs(u00) < tol
        # Note that u10 can't be zero, since u is unitary (and u00 == 0)
        gamma = angle(-u01 / u10)
        delta = angle(u01 * exp(-1im * gamma / 2))
        return 0., pi, -gamma, delta
    end
    # Handle special case if the entry (0,1) of the unitary is equal to zero
    if abs(u01) < tol
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

"""
theta, phi, lam, phase - 0.5 * (phi + lam)

- Z(phi) Y(theta) Z(lambda)
- e^{i gamma}{2})} U_3(theta,phi,lambda)
"""
@memoize function zyz_to_u3(a,b,c,d)

    # Convert exp(id)*Rz(c).Ry(b).Rz(a) angles to u3
    return (b, c, a, d -0.5*(a+c) )
end

"""
theta, phi, lam, phase - 0.5 * (phi + lam)

- Z(phi) Y(theta) Z(lambda)
- e^{i gamma}{2})} U_3(theta,phi,lambda)
"""
@memoize function u3_to_zyz(a,b,c,d)
    # Convert exp(id)*Rz(c).Ry(b).Rz(a) angles to u3
    return (b, c, a, d -0.5*(a+c) )
end

@memoize function euler_to_mat(a::Number,b::Number,c::Number,d::Number)
    """Creates unitary from angles defined by mat_to_euler function as u=exp(id)*Rz(c).Ry(b).Rz(a).
    """
    u00 = exp(-1im*(c+a)/2)*cos(b/2)
    u01 = -exp(1im*(a-c)/2)*sin(b/2)
    u10 = exp(1im*(c-a)/2)*sin(b/2)
    u11 = exp(1im*(a+c)/2)*cos(b/2)
    return [u00 u01; u10 u11].*exp(1im*d)
end

@memoize function u3_to_mat(θ::Number,ϕ::Number,λ::Number,ph::Number)
    """Creates unitary from angles defined by mat_to_euler function as u=exp(id)*Rz(c).Ry(b).Rz(a).
    """
    u00 = cos(θ/2)
    u01 = -exp(1im*λ)*sin(θ/2)
    u10 = exp(1im*ϕ)*sin(θ/2)
    u11 = exp(1im*(ϕ+λ)/2)*cos(θ/2)
    return [u00 u01; u10 u11].*exp(1im*ph)
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

function Base.:String(gc::AGateCall)
    return String(gc.gate_label);
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #


"""
    create_gates_nonparam(gate_call::AGateCall)

    Create QASM syntax for gates cached in g_map_lg dictionary. 
    Single and controlled versions generated.
    Maps X^{1/2^n} gates to H Z^{1/2^n} H

"""
function create_gates_nonparam(gate_call::AGateCall)
    if isdefined(gate_call, :ctrl)
        ct = "ctrl,tgt";
        u_pre = "c"
    else
        ct = "tgt";
        u_pre = ""
    end
    label = gate_call.gate_label
    gate = g_map_lg[label]

    # Problematic if z-gates not populated before x-gates; init to sufficient depth to prevent issues
    if String(label.gate_label)[end] == 'x'
        gate_label_sub = Symbol(String(label.gate_label)[1:end-1] * "z")
        gate = g_map[gate_label_sub]
    end

    euler_vals = GateOps.mat_to_euler(gate)
    u3_vals = GateOps.zyz_to_u3(euler_vals...)

    # QASM and X^{1/2^n} gates do not work well, so map to H Z^{1/2^n} H
    if String(label.gate_label)[end] == 'x'
        g = """gate $(u_pre)$(label.gate_label) tgt{\nh tgt;\n$(u_pre)u3($(u3_vals[1]),$(u3_vals[2]),$(u3_vals[3])) tgt;\nh tgt;\n}"""
    elseif String(label.gate_label)[end] == 'z'
        g = """gate $(u_pre)$(label.gate_label) tgt{\n$(u_pre)u3($(u3_vals[1]),$(u3_vals[2]),$(u3_vals[3])) tgt;\n}"""
    else
        error("Unsupported gate: Pauli-X and Pauli-Z currently only supported")
    end

    return g
end


"""
    create_gates_nonparam(gate_label::Symbol, num_qubits::Int
    
    Create QASM syntax for gates cached in g_map_lg dictionary. 
    Single and controlled versions generated.
    Maps X^{1/2^n} gates to H Z^{1/2^n} H

"""
function create_gates_nonparam(gate_label::Symbol, num_qubits::Int)
    ct=""
    gate = g_map[gate_label]

    if String(gate_label)[end] == 'x'
        gate_label_sub = Symbol(String(gate_label)[1:end-1] * "z")
        gate = g_map[gate_label_sub]
    end

    euler_vals = GateOps.mat_to_euler(gate)
    u3_vals = GateOps.zyz_to_u3(euler_vals...)

    # QASM and X^{1/2^n} gates do not work well, so map to H Z^{1/2^n} H
    if num_qubits == 1 && String(gate_label)[end] == 'x'
        ct = "tgt"
        g = """gate $(gate_label) tgt{\nh tgt;\nu3($(u3_vals[1]),$(u3_vals[2]),$(u3_vals[3])) tgt;\nh tgt;\n}"""
    elseif num_qubits == 2 && String(gate_label)[end] == 'x'
        ct = "ctrl,tgt"
        g = """gate c$(gate_label) $(ct){\nh tgt;\ncu3($(u3_vals[1]),$(u3_vals[2]),$(u3_vals[3])) $(ct);\nh tgt;\n}"""
    elseif num_qubits == 1 && String(gate_label)[end] == 'z'
        ct = "tgt"
        g = """gate $(gate_label) tgt{\nu3($(u3_vals[1]),$(u3_vals[2]),$(u3_vals[3])) tgt;\n}"""
    else
        ct = "ctrl,tgt"
        g = """gate c$(gate_label) $(ct){\ncu3($(u3_vals[1]),$(u3_vals[2]),$(u3_vals[3])) $(ct);\n}"""
    end

    return g
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

# QASM interop functions
function mat_to_u3(unitary::Matrix{<:Number})
    # Direct port from Qiskit https://qiskit.org/documentation/_modules/qiskit/quantum_info/synthesis/one_qubit_decompose.html

    coeff = 1/sqrt(det(unitary))
    phase = -angle(coeff)
    su_mat = coeff * unitary  # U in SU(2)

    theta = 2 * atan(abs(su_mat[2, 1]), abs(su_mat[1, 1]))
    phiplambda = 2 * angle(su_mat[2, 2])
    phimlambda = 2 * angle(su_mat[2, 1])
    phi = (phiplambda + phimlambda) / 2.0
    lam = (phiplambda - phimlambda) / 2.0
    
    return theta, phi, lam, phase - 0.5 * (phi + lam)
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

end
