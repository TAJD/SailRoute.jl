using sail_route, Test

@testset "Min angle" begin
    @test sail_route.min_angle(10.0, 20.0) ≈ 10.0
    @test sail_route.min_angle(-10.0, 30.0) ≈ 40.0
    @test sail_route.min_angle(-30.0, 10.0) ≈ 40.0
end

@testset "Resolve vectors" begin
    ws = sail_route.w_c(0.0, 10.0, 0.0, 0.0) 
    @test ws[1] ≈ 10.0
    @test ws[2] ≈ 0.0
end
 
