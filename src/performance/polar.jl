using CSV, Interpolations, DataFrames, JuMP, Ipopt

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
awa(twa, v_s, v_t) = atan(sind(twa)/(cosd(twa) + v_s/v_t))

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
    etp0 = extrapolate(itp, Flat())
    return etp0
end


"""
functionterp(itp, twa, tws)

Return interpolated performance. Convert from ms to knots here.
"""
@fastmath function perf_interp(performance::Performance, twa, tws)
    if twa < performance.polar.itp.knots[1][1]
        return 0.0
    else
        v = performance.polar(twa, tws)*performance.uncertainty
        if v > 0.0
            return v
        else
            return 0
        end
    end 
end

"""Calculate the minimum distance between two angles."""
@inline @fastmath function min_angle(a, b)
    return abs(atand(sind(a-b), cosd(a-b)))
end


"""Generic Aerrtsen typical speed loss values. Needs to be made more specific."""
function typical_aerrtsen()
    itp = interpolate(([0.2, 0.6, 1.5, 2.3, 4.2, 8.2],), [100.0, 99.0, 97.0, 94.0, 83.0, 60.0], Gridded(Linear()))
    etp = extrapolate(itp, Line())
    return etp
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
    return perf_interp(perf, awa, w_c_s)
end


"""Calculate the speed of the sailing craft given the current."""
function solve_speed_given_current(tws, twa, cs, ca, bearing, perf)
    p(ϕ) = cost_func(ϕ, tws, twa, cs, ca, perf)
    ca = wwd_to_md(ca)
    h_comp(ϕ) = p(ϕ)*sind(bearing-ϕ)-cs*sind(ca-bearing)
    v_comp(ϕ) = p(ϕ)*cosd(bearing-ϕ)
    model = Model(with_optimizer(Ipopt.Optimizer,print_level=0, warm_start_init_point="yes",max_iter=30, acceptable_tol=0.2))
    # variables
    @variable(model, 0.0 <= ϕ <= 360.0, start=bearing)         
    # register functions
    register(model, :p, 1, p, autodiff=true)
    register(model, :h_comp, 1, h_comp, autodiff=true)
    register(model, :v_comp, 1, v_comp, autodiff=true)                             
    # constraints
    @NLconstraint(model, con1, h_comp(ϕ)==0.0)                                         
    @NLconstraint(model, con2, v_comp(ϕ)>=0.0)
    # objective function
    @NLobjective(model, Max, v_comp(ϕ))
    JuMP.optimize!(model)
    return v_comp(value.(ϕ))
end


"""
Cost function to use within routing simulation.
"""
function cost_function(performance, cudi, cusp, widi, wisp, wadi, wahi, bearing)
    return solve_speed_given_current(wisp, widi, cusp, cudi, bearing, performance)
end


"""Generate range of modified polars for performance uncertainty simulation."""
function generate_performance_uncertainty_samples(polar, params, wave_m)
    unc_perf = [Performance(polar, i, 1.0, wave_m) for i in params]
end


