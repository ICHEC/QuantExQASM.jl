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

abstract type AGate end
abstract type AGateCall <: AGate end
abstract type AGateLabel <: AGate end

struct GateLabelP{IType<:Integer, NType<:Number} <: AGateLabel
    label::Symbol
    params::Union{Nothing, Dict{String, Union{IType, NType, Bool}}}
    GateLabelP{IType,NType}(gate_label) where {IType<:Integer, NType<:Number} = new(gate_label, nothing)
    GateLabelP{IType,NType}(gate_label, params) where {IType<:Integer, NType<:Number} = new(gate_label, params)
end

struct GateCall1P{IType<:Integer, NType<:Number} <: AGateCall
    gate_label::GateLabelP{IType, NType}
    target::IType
    reg::Union{String, Nothing}
    GateCall1P{IType,NType}(gate_label, target, reg) where {IType<:Integer, NType<:Number} = new(gate_label, target, reg)
    GateCall1P{IType,NType}(gate_label, target) where {IType<:Integer, NType<:Number} = new(gate_label, target, "")
end

struct GateCall2P{IType<:Integer, NType<:Number} <: AGateCall
    gate_label::GateLabelP{IType, NType}
    ctrl::IType
    target::IType
    base_gate::Union{GateCall1P{IType, NType},GateCall2P{IType, NType}, Nothing}
    reg::Union{String, Nothing}
    GateCall2P{IType,NType}(gate_label, ctrl, target, base_gate, reg) where {IType<:Integer, NType<:Number} = new(gate_label, ctrl, target, base_gate{IType,NType},reg)
    GateCall2P{IType,NType}(gate_label, ctrl, target, base_gate) where {IType<:Integer, NType<:Number} = new(gate_label, ctrl, target, base_gate{IType,NType}, "")
    GateCall2P{IType,NType}(gate_label, ctrl, target) where {IType<:Integer, NType<:Number} = new(gate_label, ctrl, target, nothing, "")
end

struct GateCallNP{IType<:Integer, NType<:Number} <: AGateCall
    gate_label::GateLabelP{IType, NType}
    ctrl::Vector{IType}
    target::IType
    base_gate::GateCall1P{IType, NType}
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

function pauli_x(q_target::Int, register::Union{String, Nothing}=nothing)
    if register == nothing
        return GateCall1(GateLabel(:x), q_target)
    else
        return GateCall1(GateLabel(:x), q_target, register)
    end
end
function pauli_y(q_target::Int, register::Union{String, Nothing}=nothing)
    if register == nothing
        return GateCall1(GateLabel(:y), q_target)
    else
        return GateCall1(GateLabel(:y), q_target, register)
    end
end
function pauli_z(q_target::Int, register::Union{String, Nothing}=nothing)
    if register == nothing
        return GateCall1(GateLabel(:z), q_target)
    else
        return GateCall1(GateLabel(:z), q_target, register)
    end
end
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

# Arbitrary U matrix can be defined with generating parameters
# Gate schema can be used to indicate construction params 
# ("r" uses given gate, sqrt depth and adjoint; "e" uses euler angles and phase)
function u(label::String, q_target::Int, params::Dict)
    return GateCall1(GateLabel(Symbol(label), params), q_target)
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

function r_x(q_target::Int, theta::Real, register::Union{String, Nothing}=nothing)
    if register == nothing
        return GateCall1(GateLabel(Symbol("r_x" * "_" * "angle=" * string(theta)), Dict("angle"=>theta)), q_target )
    else
        return GateCall1(GateLabel(Symbol("r_x" * "_" * "angle=" * string(theta)), Dict("angle"=>theta)), q_target, register )
    end
end
function r_y(q_target::Int, theta::Real, register::Union{String, Nothing}=nothing)
    if register == nothing
        return GateCall1(GateLabel(Symbol("r_y" * "_" * "angle=" * string(theta)), Dict("angle"=>theta)), q_target )
    else
        return GateCall1(GateLabel(Symbol("r_y" * "_" * "angle=" * string(theta)), Dict("angle"=>theta)), q_target, register )
    end
end
function r_z(q_target::Int, theta::Real, register::Union{String, Nothing}=nothing)
    if register == nothing
        return GateCall1(GateLabel(Symbol("r_z" * "_" * "angle=" * string(theta)), Dict("angle"=>theta)), q_target )
    else
        return GateCall1(GateLabel(Symbol("r_z" * "_" * "angle=" * string(theta)), Dict("angle"=>theta)), q_target, register )
    end
end
function r_phase(q_target::Int, theta::Real, register::Union{String, Nothing}=nothing)
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
function c_pauli_x(q_target::Int, q_ctrl::Int, register::Union{String, Nothing}=nothing)
    return GateCall2( GateLabel(:c_x), q_target, q_ctrl, GateLabel(:x), register)
end
function c_pauli_y(q_target::Int, q_ctrl::Int, register::Union{String, Nothing}=nothing)
    return GateCall2( GateLabel(:c_y), q_target, q_ctrl, GateLabel(:y), register)
end
function c_pauli_z(q_target::Int, q_ctrl::Int, register::Union{String, Nothing}=nothing)
    return GateCall2( GateLabel(:c_z), q_target, q_ctrl, GateLabel(:z), register)
end
function c_u(label::String, q_target::Int, q_ctrl::Int, params::Dict, register::Union{String, Nothing}=nothing)
    return GateCall2( GateLabel( Symbol(label), params), q_target, q_ctrl, register)
end
function c_u(label::GateLabel, q_target::Int, q_ctrl::Int, register::Union{String, Nothing}=nothing)
    return GateCall2( label, q_target, q_ctrl, register)
end

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

function c_r_x(q_target::Int, q_ctrl::Int, theta::Real, register::Union{String, Nothing}=nothing)
    return GateCall2( GateLabel( Symbol("c_r_x" * "_" * "angle=" * string(theta)), Dict("angle"=>theta) ), q_target, q_ctrl )
end
function c_r_y(q_target::Int, q_ctrl::Int, theta::Real, register::Union{String, Nothing}=nothing)
    return GateCall2( GateLabel( Symbol("c_r_y" * "_" * "angle=" * string(theta)), Dict("angle"=>theta) ), q_target, q_ctrl, Dict("theta"=>theta) )
end
function c_r_z(q_target::Int, q_ctrl::Int, theta::Real, register::Union{String, Nothing}=nothing)
    return GateCall2( GateLabel( Symbol("c_r_z" * "_" * "angle=" * string(theta)), Dict("angle"=>theta) ), q_target, q_ctrl, Dict("theta"=>theta) )
end
function c_r_phase(q_target::Int, q_ctrl::Int, theta::Real, register::Union{String, Nothing}=nothing)
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


@memoize function mat_to_euler(unitary::Matrix{<:Number})
    """Finds rotation angles (a,b,c,d) in the decomposition u=exp(id)*Rz(c).Ry(b).Rz(a).
    Direct port to Julia from qiskit: https://qiskit.org/documentation/_modules/qiskit/extensions/quantum_initializer/squ.html
    """
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


@memoize function zyz_to_u3(a,b,c,d)
    """
    theta, phi, lam, phase - 0.5 * (phi + lam)

    - Z(phi) Y(theta) Z(lambda)
    - e^{i gamma}{2})} U_3(theta,phi,lambda)
    """
    # Convert exp(id)*Rz(c).Ry(b).Rz(a) angles to u3
    return (b, c, a, d -0.5*(a+c) )
end
@memoize function u3_to_zyz(a,b,c,d)
    """
    theta, phi, lam, phase - 0.5 * (phi + lam)

    - Z(phi) Y(theta) Z(lambda)
    - e^{i gamma}{2})} U_3(theta,phi,lambda)
    """
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
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

end
