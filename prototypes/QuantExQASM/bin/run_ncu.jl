using Pkg
dir_path = @__DIR__
Pkg.activate(dir_path * "/..")
using QuantExQASM

X = QuantExQASM.NCU.GateOps.default_gates["X"]
ctrl = collect(0:2)
aux = Array{Int64,1}()
tgt = length(ctrl)

cct = QuantExQASM.NCU.apply_ncu("q", ctrl, aux, tgt, X, 0)
println(QuantExQASM.Utils.format_string_nl(cct))
