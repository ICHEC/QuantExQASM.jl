module CircuitBuilder

using .GateOps
using .GateRegistry
using DataStructures

MLLG = MutableLinkedList{<:GateOps.AGate}
gate_set = Set{GateOps.GateLabel}()

mutable struct Circuit
    ops::MLLG
    num_qubits::Int
    gate_reg::GateRegistry.Registry
end

function apply_gate(l::GateOps.GateLabel, params::Tuple)
    
end

function Base:*(cct1::Circuit, cc2::Circuit)
    append!(cct1.ops, cct2.ops)
end
function Base:+(cct1::Circuit, cc2::Circuit)
    Base:*(cct1::Circuit, cc2::Circuit)
end

function build_circuit(gate_call_list::MLLG, target_output::Symbol) 

end

function get_gate(label::String)

end

function create_gate(label::String)
    
end

end