using QuantExQASM

using Compat.Test
using TestSetExtensions
using PyCall

#include("NCU_tests.jl")

@testset "All the tests" begin
    @includetests ARGS
end

