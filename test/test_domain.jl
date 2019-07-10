using sail_route, Test

@testset "Distances" begin
    dist, bearing = sail_route.haversine(-99.436554, 41.507483, -98.315949, 38.504048)
    @test dist ≈ 187.595 atol=0.01
    dist, bearing = sail_route.euclidean(0.0, 0.0, 10.0, 10.0)
    @test dist ≈ 14.142 atol=0.01
    @test bearing ≈ 45.0 
end
