#!/usr/bin/env julia

using SailRoute, Test
@testset "Discretization error functions" begin
    value_1 = 0.970500
    value_2 = 0.968540
    value_3 = 0.961780
    h1 = 1.0
    h2 = 2.0
    h3 = 4.0
    ratio_21 = h2/h1
    ratio_32 = h3/h2
    ooc = SailRoute.order_of_convergence(value_1, value_2, value_3, ratio_21, ratio_32)
    @test ooc ≈ 1.786170 atol=4
    f_exact_21 = SailRoute.richardson_extrapolate(value_1, value_2, ratio_21, ooc)
    f_exact_32 = SailRoute.richardson_extrapolate(value_2, value_3, ratio_32, ooc)
    e21a, e21ext = SailRoute.error_estimates(value_1, value_2, f_exact_21)
    @test e21a ≈ 0.002020 atol = 5
    @test e21ext ≈ 0.000824 atol = 5
    e32a, e32ext = SailRoute.error_estimates(value_2, value_3, f_exact_32)
    gci_fine_21, gci_coarse_21 = SailRoute.gci(ratio_21, e21a, ooc)
    gci_fine_32, gci_coarse_32 = SailRoute.gci(ratio_32, e32a, ooc)
    ratio = SailRoute.asymptotic_ratio(gci_fine_21, gci_fine_32, ratio_21, ooc)
    @test ratio ≈ 0.997980 atol = 5
    # @test SailRoute.perf_interp(test_performance, 60.0, 10.0) ≈ 5.0 
    # @test SailRoute.perf_interp(test_performance, 60.0, 15.0) ≈ 7.5
    # @test SailRoute.perf_interp(test_performance, 30.0, 10.0) ≈ 0.0
end