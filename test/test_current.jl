using SailRoute, Test


function circular_performance()
    twa = [a for a in LinRange(40.0, 180.0, 10)]
    tws = [10.0, 15.0]
    perf = fill(5.0, 10, 2)
    perf[:, 2] .= 7.5
    polar = SailRoute.setup_perf_interpolation(tws, twa, perf)
    performance = SailRoute.Performance(polar, 1.0, 1.0, nothing)
    return performance
end


@testset "Performance interpolation" begin
    test_performance = circular_performance()
    @test SailRoute.perf_interp(test_performance, 60.0, 10.0) ≈ 5.0 
    @test SailRoute.perf_interp(test_performance, 60.0, 15.0) ≈ 7.5
    @test SailRoute.perf_interp(test_performance, 30.0, 10.0) ≈ 0.0
end


@testset "Cost function" begin
    @test SailRoute.cost_func(90.0, 10.0, 0.0, 0.0, 0.0, circular_performance()) ≈ 5.0
    @test SailRoute.cost_func(90.0, 10.0, 0.0, 1.0, 0.0, circular_performance()) ≈ 5.5
end


@testset "Speed given current solution" begin
    @test SailRoute.solve_speed_given_current(10.0, 0.0, 0.0, 0.0, 90.0, circular_performance()) ≈ 5.0
    @test SailRoute.solve_speed_given_current(10.0, 0.0, 1.0, 180.0, 90.0, circular_performance()) ≈ 5.408 atol=0.01
end
