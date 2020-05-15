"""
    QuantExQASM.jl - Initial implementation of OpenQASM generator for 
    different quantum algorithms.

# Example: QFT
We create a create register, labelled as `myReg`, and a range of qubit indices.
The `gen_qft` function will output a semi-colon delimited string of instructions
in OpenqASM format.

```julia
using QuantExQASM
q_register = "myReg"
qubit_indices = Array{Int64}(0:5)

gen_qft(q_register, qubit_indices)
```

To more easily visualise the output, we can format it for readability with
```julia
output = gen_qft(q_register, qubit_indices)
print(format_string_nl(output))
```
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

#@reexport using .GateOps
#include("algs/QFT.jl")
#@reexport using .QFT
#include("Utils.jl")
#@reexport using .Utils
#include("VQE.jl")
#@reexport using .VQE

#include("algs/Decomposition.jl")
#@reexport using .Decomposition
#include("algs/NCU.jl")
#@reexport using .NCU

#include("algs/Oracle.jl")
#@reexport using .Oracle
#include("algs/Diffusion.jl")
#@reexport using .Diffusion
#include("algs/Grover.jl")

end 
