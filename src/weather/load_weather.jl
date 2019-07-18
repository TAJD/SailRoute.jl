@fastmath function calc_polars(x, y)
    r = sqrt(x^2 + y^2)
    if r == 0.0
        return 0.0, 0.0
    else
        theta = atand(x, y)
        if theta < 0.0
            theta += 360.0
        end
        return r, theta
    end
end


"""Generate weather for shortest path test functions."""
function generate_constant_weather(tws_val, twa_val, cs_val, cd_val,
                          wahi_val, wadi_val, n_locs, t_len)
    empty = zeros(t_len, n_locs, n_locs)
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
