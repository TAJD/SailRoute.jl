using PyCall, Interpolations

module sail_route

include("weather/load_weather.jl")
include("performance/polar.jl")

end # module
