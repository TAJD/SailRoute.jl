#!/usr/bin/env julia


using SailRoute, Test, SafeTestsets
println("Starting tests")

@time @safetestset "Performance" begin include("test_polar.jl") end
@time @safetestset "Current" begin include("test_current.jl") end
@time @safetestset "Domain functions" begin include("test_domain.jl") end
@time @safetestset "Shortest path" begin include("test_shortest_path.jl") end
@time @safetestset "Discretization routine" begin include("test_discretization.jl") end
