using sail_route, Test


function circular_performance()
    twa = [a for a in LinRange(40.0, 180.0, 10)]
    tws = [10.0, 15.0]
    perf = fill(5.0, 10, 2)
    perf[:, 2] .= 7.5
    polar = sail_route.setup_perf_interpolation(tws, twa, perf)
    performance = sail_route.Performance(polar, 1.0, 1.0, nothing)
    return performance
end


@testset "Check performance interpolation" begin
    test_performance = circular_performance()
    @test sail_route.perf_interp(test_performance, 60.0, 10.0) ≈ 5.0 
    @test sail_route.perf_interp(test_performance, 60.0, 15.0) ≈ 7.5
    @test sail_route.perf_interp(test_performance, 30.0, 10.0) ≈ 0.0
end


@testset "Check cost function" begin
    @test sail_route.cost_func(90.0, 10.0, 0.0, 0.0, 0.0, circular_performance()) ≈ 5.0
    @test sail_route.cost_func(90.0, 10.0, 0.0, 1.0, 0.0, circular_performance()) ≈ 5.5
end

#sp = sail_route.solve_speed_given_current(10.0, 0.0, 1.0, 0.0, 90.0, circular_performance())
#println(sp)
