using CSV, Interpolations, DataFrames, Roots, Formatting, JuMP, Ipopt

"""
    awa(twa, v_s, v_t)

Calculate the apparent wind angle given true wind angle, boat speed and wind speed.

# Example
#
# ```jldoctest
# using sail_route
# c_awa = sail_route.awa(60, 3.086, 5.144)
# isapprox(0.6669807044553968, c_awa, rtol=3)
#
# # output
#
# true
# ```
# """
awa(twa::Float64, v_s::Float64, v_t::Float64) = atan(sind(twa)/(cosd(twa) + v_s/v_t))

"""
    Performance(polar, uncertainty, acceptable_failure, wave_resistance)

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
    perf = convert(Matrix, df[2:end, 2:end])
    tws = [df[1, col] for col in 2:ncol(df)]
    twa = map(x->parse(Float64, x), df[2:end, 1])
    return convert(Array{Float64}, twa), convert(Array{Float64},tws), perf
end


"""
    setup_perf_interpolation(tws, twa, perf)

Return interpolation object.
"""
function setup_perf_interpolation(tws, twa, perf)
    knots = (twa, tws)
    itp = interpolate(knots, perf, Gridded(Linear()))
    etp0 = extrapolate(itp, Line())
    #etp0 = extrapolate(itp, 0.0)
    return etp0
end


"""
functionterp(itp, twa, tws)

Return interpolated performance. Convert from ms to knots here.
"""
@inline @fastmath function perf_interp(performance::Performance, twa, tws)
    return performance.polar(twa, tws*1.94384)
end

"""Calculate the minimum distance between two angles."""
@fastmath function min_angle(a, b)
    return abs(atand(sind(a-b), cosd(a-b)))
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

function print_perf_calcs(w_c_sp, w_c_di, v_init, low_bear, high_bear)
    printfmt("w_c_sp: {1:.2f} w_c_di: {2:.2f} Bsp_init: {3:.2f}\n", w_c_sp, w_c_di, v_init)
    printfmt("low_bear: {1:.2f}, high_bear: {2:.2f}\n", low_bear, high_bear)
end


"""Generic Aerrtsen typical speed loss values. Needs to be made more specific."""
function typical_aerrtsen()
    itp = interpolate(([0.2, 0.6, 1.5, 2.3, 4.2, 8.2],), [100.0, 99.0, 97.0, 94.0, 83.0, 60.0], Gridded(Linear()))
    etp = extrapolate(itp, Line())
    return etp
end


function print_perf(awa, tws, v)
    printfmt("AWA = {1:.2f}, TWS = {2:.2f}, v = {3:.2f}\n", awa, tws, v)
end

@inline @fastmath wwd_to_md(x)=mod2pi(deg2rad(270.0-x))
@inline @fastmath md_to_wwd(x)=mod(rad2deg(x)-270.0, 360.0)

"""Resolve the wind and current vectors."""
@fastmath function w_c(twa, tws, ca, cs)
    twa = wwd_to_md(twa)
    ca = wwd_to_md(ca)
    x = tws*cos(twa) + cs*cos(ca)
    y = tws*sin(twa) + cs*sin(ca)
    w_c_s = sqrt(x^2+y^2)
    w_c_a = md_to_wwd(atan(y, x))
    return w_c_s, w_c_a
end


"""Return the speed of a sailing craft"""
@fastmath function cost_func(ϕ, tws, twa, cs, ca, perf)
    w_c_s, w_c_a = w_c(twa, tws, ca, cs)
    awa = min_angle(w_c_a, ϕ)
    bsp = perf_interp(perf, awa, w_c_s)
    return bsp
end


"""Calculate the speed of the sailing craft given the current."""
function solve_speed_given_current(tws, twa, cs, ca, bearing, perf)
    p(ϕ) = cost_func(ϕ, tws, twa, cs, ca, perf)
    h_comp(ϕ) = p(ϕ)*sind(bearing-ϕ)-cs*sind(wwd_to_md(ca)-bearing)
    v_comp(ϕ) = p(ϕ)*cos(bearing-ϕ)
    
    model = Model(with_optimizer(Ipopt.Optimizer, print_level = 0))
    
    # variables
    @variable(model, 0.0 <= ϕ <= 360.0, start=bearing)         
    # register functions
    register(model, :p, 1, p, autodiff=true)
    register(model, :h_comp, 1, h_comp, autodiff=true)
    register(model, :v_comp, 1, v_comp, autodiff=true)                             
    # constraints
    @NLconstraint(model, con, h_comp(ϕ)==0)                                         
    # objective function
    @NLobjective(model, Max, v_comp(ϕ))
    JuMP.optimize!(model)
    return p(value.(ϕ))
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
    print_perf_calcs(w_c_sp, w_c_di, v_init, resultant_0, low_bear, high_bear)
end


"""
Function to return speed without considering current.
"""
function cost_function_test(performance, cudi::Float64, cusp::Float64, widi::Float64, wisp::Float64,
                            wadi::Float64, wahi::Float64, bearing::Float64)
    return solve_speed_given_current(wisp, widi, cusp, cudi, bearing, performance)
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
   # println("W_c_di: ", w_c_di)
   # println("bearing: ", bearing)
   # println("AWA: ", min_angle(w_c_di, bearing))
    resultant(x) = hor_result(performance, w_c_di, w_c_sp, wahi, wadi, cudi, cusp, bearing, x)
    low_bearing = check_brackets(bearing-45.0, w_c_di)
    high_bearing = check_brackets(bearing+45.0, w_c_di)
    print_perf_calcs(w_c_sp, w_c_di, v, low_bearing, high_bearing)
    try
        phi = find_zero(resultant, (low_bearing, high_bearing), xatol=0.1)
        #println("phi: ", phi)
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
@fastmath cost_function(performance, cudi, cusp, widi, wisp, wadi, wahi, bearing) = cost_function_test(performance, cudi, cusp, widi, wisp, wadi, wahi, bearing)


"""Generate range of modified polars for performance uncertainty simulation."""
function generate_performance_uncertainty_samples(polar, params, wave_m)
    unc_perf = [Performance(polar, i, 1.0, wave_m) for i in params]
end

