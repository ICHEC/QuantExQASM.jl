# QuantExQASM.jl Documentation

#```@meta
#CurrentModule = QuantExQASM
#```

```@docs
QuantExQASM.Circuit.add_to_cache(label::QuantExQASM.GateOps.GateLabel, mat::Matrix{<:Number})
```

```@docs
QuantExQASM.GateOps.apply_gate_x(q_reg::String, q_idx::Union{Int, Nothing}=nothing)
QuantExQASM.GateOps.apply_gate_y(q_reg::String, q_idx::Union{Int, Nothing}=nothing)
QuantExQASM.GateOps.apply_gate_z(q_reg::String, q_idx::Union{Int, Nothing}=nothing)
QuantExQASM.GateOps.apply_gate_h(q_reg::String, q_idx::Union{Int, Nothing}=nothing)
```