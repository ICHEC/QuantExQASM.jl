module GateOps

export GateOps

struct GateCallP{IType<:Real, FType<:Real}
    gate_label::Symbol
    target::IType
    ctrl::Union{IType, Nothing}
    param::Union{FType, Nothing}
    GateCallP{IType,FType}() where {IType<:Real, FType<:Real} = new(:nothing, 0, nothing, nothing)
    GateCallP{IType,FType}(gate_label, target)  where {IType<:Real, FType<:Real} = new(gate_label, target, nothing, nothing)
    GateCallP{IType,FType}(gate_label, target, ctrl) where {IType<:Real, FType<:Real}= new(gate_label, target, ctrl, nothing)
    GateCallP{IType,FType}(gate_label, target, ctrl, param) where {IType<:Real, FType<:Real} = new(gate_label, target, ctrl, param)
end

GateCall = GateCallP{Int64, Float64}

function pauli_x(q_target::Int)
    return GateCall(:x, q_target)
end
function pauli_y(q_target::Int)
    return GateCall(:y, q_target)
end
function pauli_z(q_target::Int)
    return GateCall(:z, q_target)
end
function hadamard(q_target::Int)
    return GateCall(:h, q_target)
end

function r_x(q_target::Int, theta::Real)
    return GateCall(:r_x, q_target, nothing, theta)
end
function r_y(q_target::Int, theta::Real)
    return GateCall(:r_y, q_target, nothing, theta)
end
function r_z(q_target::Int, theta::Real)
    return GateCall(:r_z, q_target, nothing, theta)
end
function r_phase(q_target::Int, theta::Real)
    return GateCall(:r_phase, q_target, nothing, theta)
end

##
function swap(q_target::Int, q_ctrl::Int)
    return GateCall(:swap, q_target, q_ctrl)
end
function c_pauli_x(q_target::Int, q_ctrl::Int)
    return GateCall(:c_x, q_target, q_ctrl)
end
function c_pauli_y(q_target::Int, q_ctrl::Int)
    return GateCall(:c_y, q_target, q_ctrl)
end
function c_pauli_z(q_target::Int, q_ctrl::Int)
    return GateCall(:c_z, q_target, q_ctrl)
end
function c_r_x(q_target::Int, q_ctrl::Int, theta::Real)
    return GateCall(:c_r_x, q_target, q_ctrl, theta)
end
function c_r_y(q_target::Int, q_ctrl::Int, theta::Real)
    return GateCall(:c_r_y, q_target, q_ctrl, theta)
end
function c_r_z(q_target::Int, q_ctrl::Int, theta::Real)
    return GateCall(:c_r_z, q_target, q_ctrl, theta)
end
function c_r_phase(q_target::Int, q_ctrl::Int, theta::Real)
    return GateCall(:c_r_phase, q_target, q_ctrl, theta)
end
 

end
