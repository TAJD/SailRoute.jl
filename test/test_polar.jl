using SailRoute, Test

@testset "Min angle" begin
    @test SailRoute.min_angle(10.0, 20.0) ≈ 10.0
    @test SailRoute.min_angle(-10.0, 30.0) ≈ 40.0
    @test SailRoute.min_angle(-30.0, 10.0) ≈ 40.0
end

@testset "Resolve vectors" begin
    ws = SailRoute.w_c(0.0, 10.0, 0.0, 0.0) 
    @test ws[1] ≈ 10.0
    @test ws[2] ≈ 0.0
end
