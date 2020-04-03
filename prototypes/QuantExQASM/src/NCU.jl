module NCU

import ..GateOps
import ..Decomposition

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
function apply_ncu(q_reg::String, ctrls::Vector{Int}, aux::Vector{Int}, tgt::Int, gate::GateOps.Gate, call_depth::Int)
    circuit = ""
    
    local_depth::Int = call_depth + 1;
    cOps::Int = length(ctrls);

    # Use optimisied nCX for larger controls, where applicable
    if (cOps >= 5) && (length(aux) >= cOps-2) && (call_depth == 0) && (gate.label == "X")
        circuit *= NCX_Opt(q_reg, ctrls, aux, tgt)

    #Use optimised 3CU gate calls where applicable
    elseif cOps == 3
        circuit *= NCU_3Opt(q_reg, ctrls, tgt, gate)

    #Use default expansion otherwise
    elseif cOps >= 2 && cOps !=3
        circuit *= NCU_default(q_reg, ctrls, aux, tgt, gate, local_depth)
    else
        @debug "CU( $(ctrls[1]), $(tgt) )  U = $(Decomposition.u3_to_gate(gate.angles, gate.label).mat) func=apply_ncu"
        circuit *= GateOps.apply_gate_cu(gate.angles, q_reg, ctrls[1], tgt)
    end

    return circuit
end

"""

Uses additional qubits provided by `aux` register to optimise 2-qubit call depth
"""
function NCX_Opt(q_reg::String, ctrls::Vector{Int}, aux::Vector{Int}, tgt::Int)
    cct_string = ""

    @debug "CCX( $(ctrls[end]), $(aux[ 1 + length(ctrls)-3 ]), $(tgt) ) func=NCX_Opt" 
    cct_string *= GateOps.apply_gate_ccx(q_reg, ctrls[end], aux[ 1 + length(ctrls)-3 ], tgt)

    for i in reverse(2:length(ctrls)-2)
        @debug "CCX( $(ctrls[1 + i]), $(aux[1 + (i-2)]), $(aux[1 + (i-1)]) ) func=NCX_Opt"
        cct_string *= GateOps.apply_gate_ccx(q_reg, ctrls[1 + i], aux[1 + (i-2)], aux[1 + (i-1)])
    end

    @debug "CCX( $(ctrls[1]), $(ctrls[2]), $(aux[1]) ) func=NCX_Opt"
    cct_string *= GateOps.apply_gate_ccx(q_reg, ctrls[1], ctrls[2], aux[1])

    for i in 2:length(ctrls)-2
        @debug "CCX( $(ctrls[1 + i]), $(aux[1 + (i-2)]), $(aux[1 + (i-1)]) ) func=NCX_Opt"
        cct_string *= GateOps.apply_gate_ccx(q_reg, ctrls[1 + i], aux[1 + (i-2)], aux[1 + (i-1)])
    end

    @debug "CCX( $(ctrls[end]), $(aux[1 + length(ctrls) - 3]), $(tgt) ) func=NCX_Opt"
    cct_string *= GateOps.apply_gate_ccx(q_reg, ctrls[end], aux[1 + length(ctrls) - 3], tgt)

    for i in reverse(2:length(ctrls)-2)
        @debug "CCX( $(ctrls[1 + i]), $(aux[1 + (i-2)]), $(aux[1 + (i-1)]) ) func=NCX_Opt"
        cct_string *= GateOps.apply_gate_ccx(q_reg, ctrls[1 + i], aux[1 + (i-2)], aux[1 + (i-1)])
    end

    @debug "CCX( $(ctrls[1]), $(ctrls[2]), $(aux[1]) ) func=NCX_Opt"
    cct_string *= GateOps.apply_gate_ccx(q_reg, ctrls[1], ctrls[2], aux[1])

    for i in 2:length(ctrls)-2
        @debug "CCX( $(ctrls[1 + i]), $(aux[1 + (i-2)]), $(aux[1 + (i-1)]) ) func=NCX_Opt"
        cct_string *= GateOps.apply_gate_ccx(q_reg, ctrls[1 + i], aux[1 + (i-2)], aux[1 + (i-1)])
    end

    return cct_string

end

"""
Optimised gate-call depth for a 3CU gate; takes 17 gates -> 13
"""
function NCU_3Opt(q_reg::String, ctrls::Vector{Int}, tgt::Int, gate::GateOps.Gate)
    cct_string = ""

    #Requires (gate)^1/4 and sqrt(gate)^1/4 \dagger 
    g, g_adj = Decomposition.gate_root_adj(gate)
    g, g_adj = Decomposition.gate_root_adj(g)


    @debug "CU( $(ctrls[1]), $(tgt) )  U = $(Decomposition.u3_to_gate(g.angles, g.label).mat) func=NCU_3Opt"
    cct_string *= GateOps.apply_gate_cu(g.angles, q_reg, ctrls[1], tgt)
    @debug "CX( $(ctrls[1]), $(ctrls[2]) ) func=NCU_3Opt"
    cct_string *= GateOps.apply_gate_cx(q_reg, ctrls[1], ctrls[2])
    @debug "CU( $(ctrls[2]), $(tgt) )  U = $(Decomposition.u3_to_gate(g_adj.angles, g_adj.label).mat) func=NCU_3Opt"
    cct_string *= GateOps.apply_gate_cu(g_adj.angles, q_reg, ctrls[2], tgt)
    @debug "CX( $(ctrls[1]), $(ctrls[2]) ) func=NCU_3Opt"
    cct_string *= GateOps.apply_gate_cx(q_reg, ctrls[1], ctrls[2])

    @debug "CU( $(ctrls[2]), $(tgt) )  U = $(Decomposition.u3_to_gate(g.angles, g.label).mat) func=NCU_3Opt"
    cct_string *= GateOps.apply_gate_cu(g.angles, q_reg, ctrls[2], tgt)
    @debug "CX( $(ctrls[2]), $(ctrls[3]) ) func=NCU_3Opt"
    cct_string *= GateOps.apply_gate_cx(q_reg, ctrls[2], ctrls[3])
    @debug "CU( $(ctrls[3]), $(tgt) )  U = $(Decomposition.u3_to_gate(g_adj.angles, g_adj.label).mat) func=NCU_3Opt"
    cct_string *= GateOps.apply_gate_cu(g_adj.angles, q_reg, ctrls[3], tgt)
    @debug "CX( $(ctrls[1]), $(ctrls[3]) ) func=NCU_3Opt"
    cct_string *= GateOps.apply_gate_cx(q_reg, ctrls[1], ctrls[3])

    @debug "CU( $(ctrls[3]), $(tgt) )  U = $(Decomposition.u3_to_gate(g.angles, g.label).mat) func=NCU_3Opt"
    cct_string *= GateOps.apply_gate_cu(g.angles, q_reg, ctrls[3], tgt)
    @debug "CX( $(ctrls[2]), $(ctrls[3]) ) func=NCU_3Opt"
    cct_string *= GateOps.apply_gate_cx(q_reg, ctrls[2], ctrls[3])
    @debug "CU( $(ctrls[3]), $(tgt) )  U = $(Decomposition.u3_to_gate(g_adj.angles, g_adj.label).mat) func=NCU_3Opt"
    cct_string *= GateOps.apply_gate_cu(g_adj.angles, q_reg, ctrls[3], tgt)
    @debug "CX( $(ctrls[1]), $(ctrls[3]) ) func=NCU_3Opt"
    cct_string *= GateOps.apply_gate_cx(q_reg, ctrls[1], ctrls[3])
    @debug "CU( $(ctrls[3]), $(tgt) )  U = $(Decomposition.u3_to_gate(g.angles, g.label).mat) func=NCU_3Opt"
    cct_string *= GateOps.apply_gate_cu(g.angles, q_reg, ctrls[3], tgt)

    return cct_string
end

function NCU_default(q_reg::String, ctrls::Vector{Int}, aux::Vector{Int}, tgt::Int, gate::GateOps.Gate, local_depth::Int)
    cct_string = ""

    #Uses memoization to cache previously used values
    g, g_adj = Decomposition.gate_root_adj(gate)

    @debug "CU( $(ctrls[end]), $(tgt) )  U = $(Decomposition.u3_to_gate(g.angles, g.label).mat) func=NCU_default"
    cct_string *= GateOps.apply_gate_cu(g.angles, q_reg, ctrls[end], tgt)

    cct_string *= apply_ncu(q_reg, ctrls[1:end-1], vcat(aux,tgt), ctrls[end], GateOps.default_gates["X"], 0)

    @debug "CU( $(ctrls[end]), $(tgt) )  U = $(Decomposition.u3_to_gate(g_adj.angles, g_adj.label).mat) func=NCU_default"
    cct_string *= GateOps.apply_gate_cu(g_adj.angles, q_reg, ctrls[end], tgt)

    cct_string *= apply_ncu(q_reg, ctrls[1:end-1], vcat(aux,tgt), ctrls[end], GateOps.default_gates["X"], 0)
    cct_string *= apply_ncu(q_reg, ctrls[1:end-1], vcat(aux,ctrls[end]), tgt, g, local_depth )

    return cct_string
end

end
