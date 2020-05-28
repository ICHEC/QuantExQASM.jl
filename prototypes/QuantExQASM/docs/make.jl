using Documenter, QuantExQASM

# Ensure src dir is accessible
push!(LOAD_PATH,"../src/")

makedocs(
    modules = [QuantExQASM],
    clean = false,
    sitename = "QuantExQASM.jl",
    pages = Any[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Manual" => Any[
            "Gates" => Any[ "Gate Labels" => "gates/gatelabel.md", "Gate Calls" => "gates/gatecall.md" ],
            "Circuits" => "circuits/circuits.md",
            "Algorithms" => Any[ "NCU" => "algo/ncu.md", "Grover" => "algo/grover.md" ],
        ],
    ]
)