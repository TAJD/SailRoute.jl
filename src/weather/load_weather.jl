using PyCall, Interpolations

w = pyimport("routing_helper")

"Generate known weather conditions"
function sample_weather()
    wisp, widi, cusp, cudi, wahi, wadi = w[:sample_weather_scenario]()
    return wisp, widi, cusp, cudi, wahi, wadi
end

function load_dataset(path_nc, var)
    ds_var = w[:load_dataset](path_nc, var)
    return ds_var
end


"Return the wave height and direction and the wind speed and direction from an ERA5 weather file."
function process_era5_weather(path_nc, longs, lats)
    rg_wisp, rg_widi, rg_wh, rg_wd, rg_wp = w[:process_era5_weather](path_nc, longs, lats)
    return rg_wisp, rg_widi, rg_wh, rg_wd, rg_wp
end


function load_era5_weather(path_nc)
    wisp, widi, wh, wd, wp = w[:retrieve_era5_weather](path_nc)
    return wisp, widi, wh, wd, wp
end


"""Regrids the weather data to a grid which is approximately the same as the sailing domain."""
function regrid_data(ds, longs, lats)
    dataset = w.regrid_data(ds, longs, lats)
    return Array{Float64}(dataset[:values])
end


"""Regrids the weather data to a grid which is the same as the sailing domain."""
function regrid_domain(ds, req_lons, req_lats)
    values, lons, lats = w.return_data(ds)
    req_lons = mod.(req_lons .+ 360.0, 360.0) 
    interp_values = zeros((size(values)[1], size(req_lons)[1], size(req_lons)[2]))
    knots = (lats[end:-1:1], lons)
    for i in 1:size(values)[1]
        itp = interpolate(knots, values[i, end:-1:1, :], Gridded(Linear()), )
        eptl  = extrapolate(itp, Line())
        interp_values[i, end:-1:1, :] = eptl.(req_lats, req_lons)
    end
    return interp_values
end

function load_cluster(path_nc, longs, lats, var)
    ds = w[:load_cluster](path_nc, longs, lats, var)
    return ds
end


function load_era20_weather(path_nc)
    wisp, widi, wh, wd, wp, time_values = w[:retrieve_era20_weather](path_nc)
    time = [Dates.unix2datetime(Int64(i)) for i in time_values]
    return wisp, widi, wh, wd, wp, time
end


function load_era5_ensemble(path_nc, ens)
    wisp, widi, wh, wd, wp, time_values = w.retrieve_era5_ensemble(path_nc, ens)
    time = [Dates.unix2datetime(Int64(i)) for i in time_values]
    return wisp, widi, wh, wd, wp, time
end


function calc_polars(x, y)
    r = sqrt(x^2 + y^2)
    if r == 0.0
        return 0.0, Inf
    elseif x < 0.0 || y > 0.0
        theta = acosd(x/r) + 180.0
    elseif y > 0.0
        theta = acosd(x/r) 
    elseif y == 0.0
        theta = acosd(x/r) + 90.0
    elseif y < 0.0
        theta = -acosd(x/r) - 90.0
    end
    return r, abs(theta)
end


"""Load current data for mid Pacific."""
function load_current_data()
    path = ENV["HOME"]*"/sail_route_old/development/polynesian/current_data/"
    m_path = path*"meridional.csv"
    z_path = path*"zonal.csv"
    meridional = convert(Matrix{Float64}, CSV.read(m_path, delim=',', datarow=1))
    zonal = convert(Matrix{Float64}, CSV.read(z_path, delim=',', datarow=1))
    meridional[:, 1] *= (0.01*1.9438444924406)
    zonal[:, 1] *= (0.01*1.9438444924406)
    meridional = meridional[end:-1:1, :]
    zonal = zonal[end:-1:1, :]
    merid_interp = interpolate((meridional[:, 2],), meridional[:, 1], Gridded(Linear()))
    merid = extrapolate(merid_interp, Line())  
    zonal_interp = interpolate((zonal[:, 2],), zonal[:, 1], Gridded(Linear()))
    zonal = extrapolate(zonal_interp, Line())  
    lats = collect(LinRange(-25, 25, 80))
    merid_sp = merid.(lats)
    zonal_sp = zonal.(lats)
    r = zeros(size(lats))
    theta = zeros(size(lats))
    r = [calc_polars(merid_sp[i], zonal_sp[i])[1] for i in eachindex(lats)]
    theta = [calc_polars(merid_sp[i], zonal_sp[i])[2] for i in eachindex(lats)]
    r_interp = interpolate((lats,), r, Gridded(Linear()))
    r_final = extrapolate(r_interp, Line()) 
    theta_interp = interpolate((lats,), theta, Gridded(Linear()))
    theta_final = extrapolate(theta_interp, Line())
    return r_final, theta_final
end


function return_current_vectors(y, t_length)
    cusp = zeros((t_length, size(y)[1], size(y)[2]))
    cudi = zeros((t_length, size(y)[1], size(y)[2]))
    r, theta = load_current_data()
    for i in 1:t_length
        cusp[i, :, :] = r.(y)
        cudi[i, :, :] = theta.(y)
    end
    return cusp, cudi
end


"""Generate weather for shortest path test functions."""
function generate_constant_weather(tws_val, twa_val, cs_val, cd_val,
                          wahi_val, wadi_val, n_locs)
    empty = zeros(2000, n_locs, n_locs)
    tws = copy(empty)
    twa = copy(empty)
    cs = copy(empty)
    cd = copy(empty)
    wahi = copy(empty)
    wadi = copy(empty)
    for i in eachindex(empty)
        tws[i] = tws_val[1]
        twa[i] = twa_val[1]
        cs[i] = cs_val[1]
        cd[i] = cd_val[1]
        wahi[i] = wahi_val[1]
        wadi[i] = wadi_val[1]
    end
    return tws, twa, cs, cd, wahi, wadi
end
