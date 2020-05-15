module NCU

import ..GateOps
import ..Decomposition
using DataStructures

using ..CircuitBuilder

"""
    apply_ncu(q_reg::String, ctrls::Vector{Int}, aux::Vector{Int}, tgt::Int, gate_label::String, call_depth::Int)

Apply an n-qubit controlled gate operation on the given target qubit.

# Arguments
- `q_reg::String`: 
- `ctrls::Vector{Int}`: 
- `aux::Vector{Int}`: 
- `tgt::Int`:
- `gate::String`:
- `call_depth::Int`:
"""


function apply_ncu!(circuit::CircuitBuilder.Circuit, ctrls::Vector{Int}, aux::Vector{Int}, tgt::Int, gate::GateCall, call_depth::Int)

    if gate in circuit.GateSet
        #pull out gate
    else
        # add to GateSet
    end
    
    local_depth::Int = call_depth + 1;
    cOps::Int = length(ctrls);
    label = string(gate.gate_label.label)

    g = GateOps.get_gate(gate)

    gate

    # Use optimisied nCX for larger controls, where applicable
    if (cOps >= 5) && (length(aux) >= cOps-2) && (call_depth == 0) && (gate.label == Symbol("x"))
        push!(circuit.ops, NCX_Opt(q_reg, ctrls, aux, tgt) )

    #Use optimised 3CU gate calls where applicable
    elseif cOps == 3
        push!(circuit.ops, NCU_3Opt(ctrls, tgt, gate) )

    #Use default expansion otherwise
    elseif cOps >= 2 && cOps !=3
        push!(circuit.ops, NCU_default(ctrls, aux, tgt, gate, local_depth) )
    else
        @debug "CU( $(ctrls[1]), $(tgt) )  U=($(label), gate.gate_label.params)  func=apply_ncu"
        push!(circuit.ops, GateOps.apply_gate_cu(gate.angles, q_reg, ctrls[1], tgt) )
        GateOps.c_u( , tgt, ctrls[1], gate.gate_label.params)
        GateOps.c_u(string(gl.gate_label), 1, 2, gl.params)
    end

    return circuit
end

"""

Uses additional qubits provided by `aux` register to optimise 2-qubit call depth
"""
function NCX_Opt!(circuit::CircuitBuilder.Circuit, ctrls::Vector{Int}, aux::Vector{Int}, tgt::Int)

    @debug "CCX( $(ctrls[end]), $(aux[ 1 + length(ctrls)-3 ]), $(tgt) ) func=NCX_Opt" 
    push!(circuit.ops, GateOps.apply_gate_ccx(q_reg, ctrls[end], aux[ 1 + length(ctrls)-3 ], tgt) )

    for i in reverse(2:length(ctrls)-2)
        @debug "CCX( $(ctrls[1 + i]), $(aux[1 + (i-2)]), $(aux[1 + (i-1)]) ) func=NCX_Opt"
        push!(circuit.ops, GateOps.apply_gate_ccx(q_reg, ctrls[1 + i], aux[1 + (i-2)], aux[1 + (i-1)]) )
    end

    @debug "CCX( $(ctrls[1]), $(ctrls[2]), $(aux[1]) ) func=NCX_Opt"
    push!(circuit.ops, GateOps.apply_gate_ccx(q_reg, ctrls[1], ctrls[2], aux[1]) )

    for i in 2:length(ctrls)-2
        @debug "CCX( $(ctrls[1 + i]), $(aux[1 + (i-2)]), $(aux[1 + (i-1)]) ) func=NCX_Opt"
        push!(circuit.ops, GateOps.apply_gate_ccx(q_reg, ctrls[1 + i], aux[1 + (i-2)], aux[1 + (i-1)]) )
    end

    @debug "CCX( $(ctrls[end]), $(aux[1 + length(ctrls) - 3]), $(tgt) ) func=NCX_Opt"
    push!(circuit.ops, GateOps.apply_gate_ccx(q_reg, ctrls[end], aux[1 + length(ctrls) - 3], tgt) )

    for i in reverse(2:length(ctrls)-2)
        @debug "CCX( $(ctrls[1 + i]), $(aux[1 + (i-2)]), $(aux[1 + (i-1)]) ) func=NCX_Opt"
        push!(circuit.ops, GateOps.apply_gate_ccx(q_reg, ctrls[1 + i], aux[1 + (i-2)], aux[1 + (i-1)]) )
    end

    @debug "CCX( $(ctrls[1]), $(ctrls[2]), $(aux[1]) ) func=NCX_Opt"
    push!(circuit.ops, GateOps.apply_gate_ccx(q_reg, ctrls[1], ctrls[2], aux[1]) )

    for i in 2:length(ctrls)-2
        @debug "CCX( $(ctrls[1 + i]), $(aux[1 + (i-2)]), $(aux[1 + (i-1)]) ) func=NCX_Opt"
        push!(circuit.ops, GateOps.apply_gate_ccx(q_reg, ctrls[1 + i], aux[1 + (i-2)], aux[1 + (i-1)]) )
    end

    return circuit
end

"""
Optimised gate-call depth for a 3CU gate; takes 17 gates -> 13
"""
function NCU_3Opt!(circuit::CircuitBuilder.Circuit, ctrls::Vector{Int}, tgt::Int, gateLabel::Symbol)

    gate = GateOps.get_gate(gateLabel)

    #Requires (gate)^1/4 and sqrt(gate)^1/4 \dagger 
    g, g_adj = Decomposition.gate_root_adj(gate)
    g, g_adj = Decomposition.gate_root_adj(g)

    @debug "CU( $(ctrls[1]), $(tgt) )  U = $(Decomposition.u3_to_gate(g.angles, g.label).mat) func=NCU_3Opt"
    push!(circuit.ops, GateOps.apply_gate_cu(g.angles, q_reg, ctrls[1], tgt) )
    @debug "CX( $(ctrls[1]), $(ctrls[2]) ) func=NCU_3Opt"
    push!(circuit.ops, GateOps.apply_gate_cx(q_reg, ctrls[1], ctrls[2]) )
    @debug "CU( $(ctrls[2]), $(tgt) )  U = $(Decomposition.u3_to_gate(g_adj.angles, g_adj.label).mat) func=NCU_3Opt"
    push!(circuit.ops, GateOps.apply_gate_cu(g_adj.angles, q_reg, ctrls[2], tgt) )
    @debug "CX( $(ctrls[1]), $(ctrls[2]) ) func=NCU_3Opt"
    push!(circuit.ops, GateOps.apply_gate_cx(q_reg, ctrls[1], ctrls[2]) )

    @debug "CU( $(ctrls[2]), $(tgt) )  U = $(Decomposition.u3_to_gate(g.angles, g.label).mat) func=NCU_3Opt"
    push!(circuit.ops, GateOps.apply_gate_cu(g.angles, q_reg, ctrls[2], tgt) )
    @debug "CX( $(ctrls[2]), $(ctrls[3]) ) func=NCU_3Opt"
    push!(circuit.ops, GateOps.apply_gate_cx(q_reg, ctrls[2], ctrls[3]) )
    @debug "CU( $(ctrls[3]), $(tgt) )  U = $(Decomposition.u3_to_gate(g_adj.angles, g_adj.label).mat) func=NCU_3Opt"
    push!(circuit.ops, GateOps.apply_gate_cu(g_adj.angles, q_reg, ctrls[3], tgt) )
    @debug "CX( $(ctrls[1]), $(ctrls[3]) ) func=NCU_3Opt"
    push!(circuit.ops, GateOps.apply_gate_cx(q_reg, ctrls[1], ctrls[3]) )

    @debug "CU( $(ctrls[3]), $(tgt) )  U = $(Decomposition.u3_to_gate(g.angles, g.label).mat) func=NCU_3Opt"
    push!(circuit.ops, GateOps.apply_gate_cu(g.angles, q_reg, ctrls[3], tgt) )
    @debug "CX( $(ctrls[2]), $(ctrls[3]) ) func=NCU_3Opt"
    push!(circuit.ops, GateOps.apply_gate_cx(q_reg, ctrls[2], ctrls[3]) )
    @debug "CU( $(ctrls[3]), $(tgt) )  U = $(Decomposition.u3_to_gate(g_adj.angles, g_adj.label).mat) func=NCU_3Opt"
    push!(circuit.ops, GateOps.apply_gate_cu(g_adj.angles, q_reg, ctrls[3], tgt) )
    @debug "CX( $(ctrls[1]), $(ctrls[3]) ) func=NCU_3Opt"
    push!(circuit.ops, GateOps.apply_gate_cx(q_reg, ctrls[1], ctrls[3]) )
    @debug "CU( $(ctrls[3]), $(tgt) )  U = $(Decomposition.u3_to_gate(g.angles, g.label).mat) func=NCU_3Opt"
    push!(circuit.ops, GateOps.apply_gate_cu(g.angles, q_reg, ctrls[3], tgt) )

    return circuit
end

function NCU_default!(circuit::CircuitBuilder.Circuit, ctrls::Vector{Int}, aux::Vector{Int}, tgt::Int, gateCall::GateOps.GateCall, local_depth::Int)

    gate = GateOps.get_gate(gateLabel)
    gateLabel_sqrt = replace(gateLabel, "sqrt="=>q_reg)

    #Uses memoization to cache previously used values
    g, g_adj = Decomposition.gate_root_adj(gate)
    GateOps.add_gate!(GateOps.GateSet, gateLabel_sqrt, gate::Matrix{<:Number})

    @debug "CU( $(ctrls[end]), $(tgt) )  U = $(Decomposition.u3_to_gate(g.angles, g.label).mat) func=NCU_default"
    push!(circuit.ops, GateOps.apply_gate_cu(g.angles, q_reg, ctrls[end], tgt) )

    push!(circuit.ops, ctrls[1:end-1], vcat(aux,tgt), ctrls[end], GateOps.default_gates["X"], 0)

    @debug "CU( $(ctrls[end]), $(tgt) )  U = $(Decomposition.u3_to_gate(g_adj.angles, g_adj.label).mat) func=NCU_default"
    push!(circuit.ops, GateOps.apply_gate_cu(g_adj.angles, q_reg, ctrls[end], tgt) )

    apply_ncu!(circuit, ctrls[1:end-1], vcat(aux,tgt), ctrls[end], GateOps.default_gates["X"], 0)
    apply_ncu!(circuit, ctrls[1:end-1], vcat(aux,ctrls[end]), tgt, g, local_depth )

    return circuit
end

end
