using CSV, Interpolations, DataFrames, Roots, PyCall

machinery = pyimport("importlib.machinery")
filename = abspath(joinpath(dirname(@__FILE__),"generate_performance.py"))
loader = machinery.SourceFileLoader("performance",filename)
perf = loader.load_module("performance")

"""
    Performance(twa, tws, boat_perf)

Type to hold sailing craft performance information.

Arguments:

- twa: Array of true wind angles for each point in polar plot.
- tws: Array of true wind speeds for each point in polar plot.
- boat_perf: Interpolations instance.

"""
struct Performance
    polar
    uncertainty::Float64
    acceptable_failure::Float64
    wave_resistance
end


"""
    load_file(path)

Load the file specified by path.
"""
function load_file(path)
    df = CSV.read(path, delim=';', datarow=1)
    perf = convert(Array{Float64}, df[2:end, 2:end])
    tws = convert(Array{Float64}, df[1, 2:end])
    twa = map(x->parse(Float64, x), df[2:end, 1])
    return twa, dropdims(tws, dims=1), perf
end


"""
    generate_performance(lower_twa, x_tws, ratio)


generate_performance(120.0, [0.0, 5.0, 10.0, 20.0, 21.0], 0.3)

Generate synthetic sailing craft performance data
"""
function generate_performance(lower_twa, x_tws, ratio)
    x_tws, xi, performance = perf.generate_performance(lower_twa, x_tws, ratio)
    return xi, x_tws, performance
end

"""
    setup_perf_interpolation(tws, twa, perf)

Return interpolation object.
"""
function setup_perf_interpolation(tws, twa, perf)
    knots = (twa, tws)
    itp = interpolate(knots, perf, Gridded(Linear()))
    # etp0 = extrapolate(itp, Line())
    etp0 = extrapolate(itp, 0.0)
    return etp0
end



"""Generic Aerrtsen typical speed loss values. Needs to be made more specific."""
function typical_aerrtsen()
    itp = interpolate(([0.2, 0.6, 1.5, 2.3, 4.2, 8.2],), [100.0, 99.0, 97.0, 94.0, 83.0, 60.0], Gridded(Linear()))
    etp = extrapolate(itp, Line())
    return etp
end

"""
    perf_interp(itp, twa, tws)

Return interpolated performance. Convert from ms to knots here.
"""
@inline @fastmath function perf_interp(performance::Performance, twa, tws, wahi, wadi)
    # if performance.wave_resistance == nothing
    #     @show "Wave resistance" wave_res = 1.0
    # else
    #     @show wave_res = performance.wave_resistance(wahi)/100.0
    # end
    wave_res = 1.0
    return performance.polar(twa, tws*1.94384)*performance.uncertainty*wave_res
end


"""Calculate the horizontal component of the current."""
@inline @fastmath function h(cudi, cusp, bearing)
    cusp*sind(cudi-bearing)
end


"""Calculate the resultant of the horizontal components of the boat and the current. The horizontal components of the speed of the vessel due to heading change and the current must be equal.

Bearing is the angle between waypoints. Heading is the actual course followed.

"""
@inline @fastmath function hor_result(performance, w_c_di, w_c_sp, wahi, wadi, cudi, cusp, bearing, heading)
    v = perf_interp(performance, min_angle(w_c_di, bearing), w_c_sp, wahi, wadi)
    return v*sind(bearing-heading) - cusp*sind(cudi-bearing)
end


"""Check if the awa suggested is within the possible values for awa from the polar data."""
@inline @fastmath function check_brackets(bearing, twa)
    awa = min_angle(bearing, twa[1])
    heading_awa = awa
    max_awa = 170.0
    min_awa = 20.0
    max_awa_delta = max_awa - awa
    min_awa_delta = awa - min_awa
    delta = 0.0
    if min_awa_delta < 0 && max_awa_delta > 0
        bearing -= min_awa_delta 
    elseif min_awa_delta < 0 && max_awa_delta < 0
        bearing += max_awa_delta
    end
    return mod(bearing, 360.0)
end


"""
cost_function_canoe(performance, cudi::Float64, cusp::Float64,
                       widi::Float64, wisp::Float64,
                       wadi::Float64, wahi::Float64,
                       bearing::Float64)


Calculate the speed of the sailing craft given the failure model and environmental conditions. Don't forget the rename function below! 
"""
function cost_function_canoe(performance::Performance,
                       cudi, cusp,
                       widi, wisp,
                       wadi, wahi,
                       bearing)
    w_c_di = mod(widi + cudi, 360.0)
    w_c_sp = wisp + cusp
    v = perf_interp(performance, min_angle(w_c_di, bearing), w_c_sp, wahi, wadi)
    resultant(x) = hor_result(performance, w_c_di, w_c_sp, wahi, wadi, cudi, cusp, bearing, x)
    low_bearing = check_brackets(bearing-45.0, w_c_di)
    high_bearing = check_brackets(bearing+45.0, w_c_di)
    try 
        phi = find_zero(resultant, (low_bearing, high_bearing), xatol=0.1)
        v = perf_interp(performance, min_angle(w_c_di, phi), w_c_sp, wahi, wadi)
        if v + cusp < 0.0  # if the speed is less than 0 then it is going backwards - return 0 
            return 0.0
        elseif (min_angle(bearing, wadi) < 40.0) && (wahi > 0.2)
            return 0.5*(v+cusp) # reduce the speed by half if heading into waves
        else
            return v + cusp
        end
    catch ArgumentError
        return 0.0
    end
end


"""
cost_function_conventional(performance, cudi::Float64, cusp::Float64,
                       widi::Float64, wisp::Float64,
                       wadi::Float64, wahi::Float64,
                       bearing::Float64)


Calculate the speed of the sailing craft given the failure model and environmental conditions. Don't forget the rename function below! 
"""
function cost_function_conventional(performance::Performance,
                       cudi, cusp,
                       widi, wisp,
                       wadi, wahi,
                       bearing)
    w_c_di = mod(widi + cudi, 360.0)
    w_c_sp = wisp + cusp
    v = perf_interp(performance, min_angle(w_c_di, bearing), w_c_sp, wahi, wadi)
    resultant(x) = hor_result(performance, w_c_di, w_c_sp, wahi, wadi, cudi, cusp, bearing, x)
    low_bearing = check_brackets(bearing-45.0, w_c_di)
    high_bearing = check_brackets(bearing+45.0, w_c_di)
    try 
        phi = find_zero(resultant, (low_bearing, high_bearing), xatol=0.1)
        v = perf_interp(performance, min_angle(w_c_di, phi), w_c_sp, wahi, wadi)
        if min_angle(phi, bearing) > 60.0
            return 0.0
        end
        if v + cusp < 0.0
            return 0.0
        else
            return v + cusp
        end
    catch ArgumentError
        return 0.0
    end
end

# this wrapper determines the performance function used in the simulation
@fastmath cost_function(performance, cudi, cusp, widi, wisp, wadi, wahi, bearing) = cost_function_canoe(performance, cudi, cusp, widi, wisp, wadi, wahi, bearing)


"""Generate range of modified polars for performance uncertainty simulation."""
function generate_performance_uncertainty_samples(polar, params, wave_m)
    unc_perf = [Performance(polar, i, 1.0, wave_m) for i in params]
end
