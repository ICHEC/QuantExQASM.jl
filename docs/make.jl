using Documenter, QuantExQASM

makedocs(
    modules = [QuantExQASM],
    clean = false,
    sitename = "QuantExQASM.jl",
    pages = Any[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Manual" => Any[
            "Gates" => "gates_circuits/GateOps.md",
            "Circuits" => "gates_circuits/circuits.md",
            "Algorithms" => Any[ "NCU" => "algo/ncu.md", "Grover" => "algo/grover.md" ],
        ],
    ]
)
deploydocs(
    repo = "github.com/ICHEC/QuantExQASM.jl.git",
)
