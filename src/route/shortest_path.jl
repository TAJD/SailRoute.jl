using Dates

struct Route
    lon1::Float64
    lon2::Float64
    lat1::Float64
    lat2::Float64
    x_nodes::Int64
    y_nodes::Int64
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
        append!(X, x_locs[idx])
        append!(Y, y_locs[idx])
    end
    return hcat(X, Y) 
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


"""Convert hours in float to an index value to be added to the start index. This needs to be changed based on the weather data used."""
@inline function convert_time(time::Float64)
    mm = round(time)
    return Dates.Hour(mm)
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
        d, b = haversine(route.lon1, route.lat1, x[1, idx], y[1, idx])
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
        if speed >= 0.0
            earliest_times[1, idx] = d/speed
        end
    end
    for idy in 1:idy_range-1
        for idx1 in 1:idx_range
            if isinf(earliest_times[idy, idx1]) == false
                t = start_time + convert_time(earliest_times[idy, idx1])
                t_idx = time_to_index(t, times)
                wd_int = widi[t_idx, idx1, idy]
                ws_int = wisp[t_idx, idx1, idy]
                wadi_int = wadi[t_idx, idx1, idy]
                wahi_int = wahi[t_idx, idx1, idy]
                cs_int = cusp[t_idx, idx1, idy]
                cd_int = cudi[t_idx, idx1, idy]
                @simd for idx2 in 1:idx_range
                    d, b = haversine(x[idy, idx1], y[idy, idx1],
                                        x[idy+1, idx2], y[idy+1, idx2])
                    speed = cost_function(performance, cd_int, cs_int,
                                                    wd_int, ws_int, wadi_int, wahi_int, b)
                    if speed >= 0.0
                        tt = earliest_times[idy, idx1] + d/speed
                        if earliest_times[idy+1, idx2] > tt
                            earliest_times[idy+1, idx2] = tt
                            prev_node[idy+1, idx2] = node_indices[idy, idx1]
                        end
                    end
                end
            end
        end
    end
    
    @simd for idx in 1:idx_range
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
            if speed >= 0.0 
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
    x_path = vcat([route.lon1], locs[:, 1], [route.lon2])
    y_path = vcat([route.lat1], locs[:, 2], [route.lat2])
    locs = hcat(x_path, y_path)
    return arrival_time, locs, earliest_times
end