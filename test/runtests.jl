using QuantExQASM
using PicoQuant

using Compat.Test
using TestSetExtensions
using PyCall

include("test_utils.jl")

@testset "All the tests" begin
    @includetests ARGS
end

