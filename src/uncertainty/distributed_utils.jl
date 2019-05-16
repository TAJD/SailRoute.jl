using Distributed, SharedArrays, CSV, Interpolations, Dates


"""Create a custom iterator which breaks up a range based on the processor number"""
function myrange(q::SharedArray) 
    @show idx = indexpids(q)
    if idx == 0 # This worker is not assigned a piece
        return 1:0, 1:0
    end
    nchunks = length(procs(q))
    splits = [round(Int, s) for s in range(0, stop=size(q,2), length=nchunks+1)]
    1:size(q,1), splits[idx]+1:splits[idx+1]
end


function route_solve_chunk!(results, t_range, p_range, 
                            sim_times, perfs,
                            route, time_indexes, x, y, 
                            wisp, widi, wadi, wahi, cusp, cudi)
    @show t_range, p_range
    for t in t_range, p in p_range
        output = sail_route.route_solve(route, perfs[p], sim_times[t], time_indexes, x, y, wisp, widi, wadi, wahi, cusp, cudi)
        @show results[t, p] = output[1]
        output = nothing
    end
end


function route_solve_save_path_chunk!(results, t_range, p_range, 
                                      sim_times, perfs, 
                                      x_results, y_results, et_results,
                                      route, time_indexes, x, y,
                                      wisp, widi, wadi, wahi, cusp, cudi)
    @show t_range, p_range
    for t in t_range, p in p_range
        output = sail_route.route_solve(route, perfs[p], sim_times[t],                                                                            time_indexes, x, y, wisp, widi, wadi, wahi, cusp, cudi)
        if isinf(output[1]) == true
            @show output[1]
            continue
        else 
            @show results[t, p] = output[1]
            x_results[t, p, :] = output[2][:, 1]
            y_results[t, p, :] = output[2][:, 2]
            et_results[t, p, :, :] = output[3]
            output = nothing
        end
    end
end


route_solve_shared_sp_chunk!(results, sim_times, perfs, x_results, y_results, et_results, route, time_indexes, x, y, wisp, widi, wadi, wahi, cusp, cudi) = route_solve_save_path_chunk!(results, myrange(results)..., sim_times, perfs, x_results, y_results, et_results, route, time_indexes, x, y, wisp, widi, wadi, wahi, cusp, cudi)

