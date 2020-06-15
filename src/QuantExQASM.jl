"""
    QuantExQASM.jl - implementation of circuit generator for 
    different quantum algorithms.
"""
module QuantExQASM

#Automagically export all symbols from included modules
using Reexport

include("GateOps.jl")
include("CircuitList.jl")
include("Circuit.jl")
include("algs/NCU.jl")
include("algs/Oracle.jl")
include("algs/Diffusion.jl")
include("algs/Grover.jl")

end 
