module SailRoute

export
	generate_performance,
	awa

include("uncertainty/discretization_error.jl")
include("route/domain.jl")
include("performance/polar.jl")
include("route/shortest_path.jl")
include("weather/load_weather.jl")

end # module
