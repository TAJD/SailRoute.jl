#!/usr/bin/env julia


using sail_route, Test, SafeTestsets
println("Starting tests")


@time @safetestset "Polar tests" begin include("test_polar.jl") end
