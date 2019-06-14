#!/usr/bin/env julia


using sail_route, Test, SafeTestsets
println("Starting tests")


@time @safetestset "Performance" begin include("test_polar.jl") end
@time @safetestset "Current" begin include("test_current.jl") end
@time @safetestset "Domain functions" begin include("test_domain.jl") end
