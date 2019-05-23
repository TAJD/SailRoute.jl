using Dates, PyCall, Formatting

dt = pyimport("datetime")

struct Route
    lon1::Float64
    lon2::Float64
    lat1::Float64
    lat2::Float64
    x_nodes::Int64
    y_nodes::Int64
end


"""Check if an array is monotonic. Works for both directions."""
function check_monotonic(array)                                                  
    return all(array[i] .> array[i+1] for i in range(1, length=length(array)-1))
end


"""Identify the shortest path given arrays of locations and an array of the earliest time at each point."""
function shortest_path(indx, pindx, sp)
    ix = findfirst(isequal(sp[end]), indx)
    if ix == nothing
        return sp
    end
    pix = pindx[ix]
    append!(sp, pix)
    if pix == 0.0
        return Array(sp)
    else
        return shortest_path(indx, pindx, sp)
    end
end


"""Get locations of the shortest path"""
function get_locs(indx, sp, x_locs, y_locs)
    X = []
    Y = []
    for k in sp[1:end-1]
        idx = findfirst(isequal(k), indx)
        # x_locs[idx]
        # y_locs[idx]
        append!(X, x_locs[idx])
        append!(Y, y_locs[idx])
    end
    return hcat(X, Y) 
end


"""Get the earliest times for a given path"""
function get_path_et(indx, sp, earliest_times)
    ET = []
    for k in sp[1:end-1]
        idx = findfirst(isequal(k), indx)
        append!(ET, earliest_times[idx])
    end
    return ET
end


"Find index of initial start time in array"
@inline function time_to_index(time, time_values)
    idx = findfirst(isequal(time), time_values)
    if isequal(idx, nothing) == true
        return size(time_values)[1]
    else
        return idx 
    end
end


"""Convert hours in float to an index value to be added to the start index."""
@inline function convert_time(time::Float64)
    mm = round(time*60/320)
    return Dates.Hour(mm*3.0)
end

function print_env(wisp, widi, wahi, wadi, cusp, cudi, sp)
    printfmt("Wisp {1:.2f} Widi {2:.2f} Wahi {3:.2f} Wadi {4:.2f} Cusp {5:.2f} Cudi {6:.2f} Bsp {7:.2f} \n", wisp, widi, wahi, wadi, cusp, cudi, sp)
end


"Time dependent shortest path."
function route_solve(route::Route, performance::Performance, start_time::DateTime, times, x, y, wisp, widi, wadi, wahi, cusp, cudi)
    start_time_idx = time_to_index(start_time, times)
    earliest_times = fill(Inf, size(x))
    prev_node = zero(x)
    node_indices = reshape(1:length(x), size(x))
    arrival_time = Inf
    final_node = 0
    earliest_times = fill(Inf, size(x))
    idx_range = size(x)[2]
    idy_range = size(x)[1]
    @simd for idx in 1:idx_range
        @inbounds d, b = haversine(route.lon1, route.lat1, x[1, idx], y[1, idx])
        wd_int = widi[start_time_idx, idx, 1]
        ws_int = wisp[start_time_idx, idx, 1]
        cs_int = cusp[start_time_idx, idx, 1]
        cd_int = cudi[start_time_idx, idx, 1]
        wadi_int = wadi[start_time_idx, idx, 1]
        if isnan(wadi_int) == true
            wadi_int = cd_int
        end
        wahi_int = wahi[start_time_idx, idx, 1]
        if isnan(wahi_int) == true
            wahi_int = 0.0
        end
        speed = cost_function(performance, cd_int, cs_int, wd_int, ws_int,
                              wadi_int, wahi_int, b)
        if speed != Inf
            earliest_times[1, idx] = d/speed
        end
    end
    @inbounds for idy in 1:idy_range-1
        @inbounds for idx1 in 1:idx_range
            if isinf(earliest_times[idy, idx1]) == false
                @inbounds t = start_time + convert_time(earliest_times[idy, idx1])
                t_idx = time_to_index(t, times)
                wd_int = widi[t_idx, idx1, idy]
                ws_int = wisp[t_idx, idx1, idy]
                wadi_int = wadi[t_idx, idx1, idy]
                wahi_int = wahi[t_idx, idx1, idy]
                cs_int = cusp[t_idx, idx1, idy]
                cd_int = cudi[t_idx, idx1, idy]
                @simd for idx2 in 1:idx_range
                    @inbounds d, b = haversine(x[idy, idx1], y[idy, idx1],
                                        x[idy+1, idx2], y[idy+1, idx2])
                    @inbounds speed = cost_function(performance, cd_int, cs_int,
                                                    wd_int, ws_int, wadi_int, wahi_int, b)
                    tt = earliest_times[idy, idx1] + d/speed
                    if earliest_times[idy+1, idx2] > tt
                        earliest_times[idy+1, idx2] = tt
                        prev_node[idy+1, idx2] = node_indices[idy, idx1]
                    end
                end
            end
        end
    end
    
    @inbounds @simd for idx in 1:idx_range
        if isinf(earliest_times[end, idx]) == false
            d, b = haversine(x[end, idx], y[end, idx], route.lon2, route.lat2)
            t = start_time + convert_time(earliest_times[end, idx])
            t_idx = time_to_index(t, times)
            wd_int = widi[t_idx, idx, end]
            ws_int = wisp[t_idx, idx, end]
            wadi_int = wadi[t_idx, idx, end]
            wahi_int = wahi[t_idx, idx, end]
            cs_int = cusp[t_idx, idx, end]
            cd_int = cudi[t_idx, idx, end]
            speed = cost_function(performance, cd_int, cs_int,
                                wd_int, ws_int, wadi_int, wahi_int, b)
            if speed != Inf 
                tt = earliest_times[end, idx] + d/speed
                if arrival_time > tt
                    arrival_time = tt
                    final_node = node_indices[end, idx]
                end
            end
        end
    end
    sp = shortest_path(node_indices, prev_node, [final_node])
    locs = get_locs(node_indices, sp, x, y)
    return arrival_time, locs, earliest_times
end



"""Numerical error reduction routine applied to real polynesian weather scenarios."""
function poly_discretization_routine(route, perf, start_time, times, wisp, widi, wahi, wadi)
    d, b = haversine(route.lon1, route.lat1, route.lon2, route.lat2)
    d_n_range = [5.0, 2.5, 1.0] # normalized height
    results = []
    gci = []
    extrap = []
    ooc = []
    for i in d_n_range # Loop 1
        min_dist = d*i/100.0
        n = calc_nodes(route.lon1, route.lon2, route.lat1, route.lat2, min_dist)
        sim_route = Route(route.lon1, route.lon2, route.lat1, route.lat2, n, n)
        x, y, wisp_i, widi_i, wadi_i, wahi_i= generate_inputs(sim_route, wisp, widi,
                                                              wadi, wahi)
        dims = size(wisp_i)
        cusp_i, cudi_i = return_current_vectors(y, dims[1])
        res = route_solve(sim_route, perf,
                                     start_time, 
                                     times, 
                                     x, y,
                                     wisp_i, widi_i,
                                     wadi_i, wahi_i,
                                     cusp_i, cudi_i)
        push!(results, res[1])
    end
    gci = GCI_calc(results[end], results[end-1], results[end-2],
                   d_n_range[end], d_n_range[end-1], d_n_range[end-2])
    extrap = extrap_value(results[end], results[end-1], results[end-2],
                          d_n_range[end], d_n_range[end-1], d_n_range[end-2])
    ooc = ooc_value(results[end], results[end-1], results[end-2],
                    d_n_range[end], d_n_range[end-1], d_n_range[end-2])
    if ooc > 1.0 && check_monotonic(results[end-2:end])==true && gci < 0.05
        return extrap, gci, ooc
    else 
        iterations = 0
        while ooc < 1.0  || check_monotonic(results[end-2:end])==false || gci > 0.05
            if float(iterations) > 2.0
                return extrap, gci
            end
            iterations += 1
            d_n_m = 3.0*d_n_range[end]/4.0 # calculate next multiple of d_n
            push!(d_n_range, d_n_m)
            min_dist = d*d_n_m/100.0
            n = calc_nodes(route.lon1, route.lon2, route.lat1, route.lat2, min_dist)
            sim_route = Route(route.lon1, route.lon2, route.lat1, route.lat2, n, n)
            x, y, wisp_i, widi_i, wadi_i, wahi_i = generate_inputs(sim_route, wisp, widi, wadi, wahi)
            dims = size(wisp_i)
            cusp_i, cudi_i = return_current_vectors(y, dims[1])
            res = route_solve(sim_route, perf,
                                        start_time, 
                                        times, 
                                        x, y,
                                        wisp_i, widi_i,
                                        wadi_i, wahi_i,
                                        cusp_i, cudi_i)
            push!(results, res[1])
            gci = GCI_calc(results[end], results[end-1], results[end-2],
                d_n_range[end], d_n_range[end-1], d_n_range[end-2])
            extrap = extrap_value(results[end], results[end-1], results[end-2],
                        d_n_range[end], d_n_range[end-1], d_n_range[end-2])
            ooc = ooc_value(results[end], results[end-1], results[end-2],
                    d_n_range[end], d_n_range[end-1], d_n_range[end-2])
        end
        return extrap, gci, ooc
    end
end


"""Routing over a flat surface."""
function cartesian_route_solve(route::Route, performance::Performance, start_time::DateTime,
                               times, x, y,
                               wisp, widi,
                               wadi, wahi,
                               cusp, cudi)
    start_time_idx = time_to_index(start_time, times)
    earliest_times = fill(Inf, size(x))
    prev_node = zero(x)
    node_indices = reshape(1:length(x), size(x))
    arrival_time = Inf
    final_node = 0
    earliest_times = fill(Inf, size(x))
    for idx in 1:size(x)[2]
        d, b = euclidean(route.lon1, route.lat1, x[1, idx], y[1, idx])
        wd_int = widi[start_time_idx, idx, 1]
        ws_int = wisp[start_time_idx, idx, 1]
        wadi_int = wadi[start_time_idx, idx, 1]
        wahi_int = wahi[start_time_idx, idx, 1]
        cs_int = cusp[start_time_idx, idx, 1]
        cd_int = cudi[start_time_idx, idx, 1]
        speed = cost_function(performance, cd_int, cs_int, wd_int, ws_int,
                    wadi_int, wahi_int, b)
        if speed != Inf
            earliest_times[1, idx] = d/speed
        else
            earliest_times[1, idx] = Inf
        end
    end
    @inbounds for idy in 1:size(x)[1]-1
        @inbounds for idx1 in 1:size(x)[2]
            if isinf(earliest_times[idy, idx1]) == false
                @inbounds t = start_time + convert_time(earliest_times[idy, idx1])
                t_idx = time_to_index(t, times)
                wd_int = widi[t_idx, idx1, idy]
                ws_int = wisp[t_idx, idx1, idy]
                wadi_int = wadi[t_idx, idx1, idy]
                wahi_int = wahi[t_idx, idx1, idy]
                cs_int = cusp[t_idx, idx1, idy]
                cd_int = cudi[t_idx, idx1, idy]
                @inbounds for idx2 in 1:size(x)[2]
                    d, b = euclidean(x[idy, idx1], y[idy, idx1],
                                        x[idy+1, idx2], y[idy+1, idx2])
                    @inbounds speed = cost_function(performance, cd_int, cs_int,
                                                    wd_int, ws_int, wadi_int, wahi_int, b)
                    tt = earliest_times[idy, idx1] + d/speed
                    if earliest_times[idy+1, idx2] > tt
                        earliest_times[idy+1, idx2] = tt
                        prev_node[idy+1, idx2] = node_indices[idy, idx1]
                    end
                    earliest_times[idy+1, idx2]
                end
            end
        end
    end

    @inbounds for idx in 1:size(x)[2]
        if isinf(earliest_times[end, idx]) == false
            d, b = euclidean(x[end, idx], y[end, idx], route.lon2, route.lat2)
            t = start_time + convert_time(earliest_times[end, idx])
            t_idx = time_to_index(t, times)
            wd_int = widi[t_idx, idx, end]
            ws_int = wisp[t_idx, idx, end]
            wadi_int = wadi[t_idx, idx, end]
            wahi_int = wahi[t_idx, idx, end]
            cs_int = cusp[t_idx, idx, end]
            cd_int = cudi[t_idx, idx, end]
            speed = cost_function(performance, cd_int, cs_int,
                                  wd_int, ws_int, wadi_int, wahi_int, b)
            if speed != Inf 
                tt = earliest_times[end, idx] + d/speed
                if arrival_time > tt
                    arrival_time = tt
                final_node = node_indices[end, idx]
                end
            end
        end
    end
    sp = shortest_path(node_indices, prev_node, [final_node])
    locs = get_locs(node_indices, sp, x, y)
    return arrival_time, locs, earliest_times, x, y
end



"""Run cartesian shortest path algorithm over the range of non dimensional grid heights given in d_n_range"""
function constant_weather_discretization_routine(route, perf, start_time, times, d_n_range)
    results = []
    d, b = haversine(route.lon1, route.lat1, route.lon2, route.lat2)
    for i in d_n_range # Loop 1
        d_n = i/100.0 * d
        n = calc_nodes(route.lon1, route.lon2, route.lat1, route.lat2, d_n)
        lats = LinRange(0.0, 50.0, n)
        lons = LinRange(0.0, 50.0, n)
        x = [i for i in lons, j in lats]
        y = [j for i in lons, j in lats]
        wisp_i, widi_i, cusp_i, cudi_i, wahi_i, wadi_i = generate_constant_weather(10.0, 0.0, 0.0, 0.0, 0.0, 0.0, n)
        sim_route = Route(route.lon1, route.lon2, route.lat1, route.lat2, n, n)
        res = cartesian_route_solve(sim_route, perf,
                                     start_time, 
                                     times, 
                                     x, y,
                                     wisp_i, widi_i,
                                     wadi_i, wahi_i,
                                     cusp_i, cudi_i)
        push!(results, res[1])
    end
    return results
end
