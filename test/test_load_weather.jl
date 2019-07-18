using SailRoute, Test


@testset "Test conversion to polar coordinates" begin
    r_1, theta_1 = SailRoute.calc_polars(1, 1)
    @test r_1 ≈ sqrt(2)
    @test theta_1 ≈ 45.0
    r_2, theta_2 = SailRoute.calc_polars(1, -1)
    @test theta_2 ≈ 135
    r_3, theta_3 = SailRoute.calc_polars(-1, -1)
    @test theta_3 ≈ 225
    r_4, theta_4 = SailRoute.calc_polars(-1, 1)
    @test theta_4 ≈ 315
end


@testset "Generate constant weather sets" begin
    tws_val = 1.0; twa_val = 90.0; cs_val=0.0; cd_val = 90.0;
    wahi_val = 0.0; wadi_val=90.0; n_locs=16;
    t_len = 100
    tws, twa, cs, cd, wahi, wadi = SailRoute.generate_constant_weather(tws_val, twa_val, cs_val, cd_val,
    wahi_val, wadi_val, n_locs, t_len)
    @test tws[1, 1, 1] ≈ tws_val
    @test twa[t_len, n_locs, 2] ≈ twa_val
    @test cs[5, 3, n_locs] ≈ cs_val
    @test cd[1, 1, 1] ≈ cd_val
    @test wahi[1, 1, 1] ≈ wahi_val
    @test wadi[1, 1, 1] ≈ wadi_val
end
