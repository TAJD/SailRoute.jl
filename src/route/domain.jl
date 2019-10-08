
"""
    haversine(lon1, lat1, lon2, lat2)

Calculate the haversine distance and bearing. Distance is in nm.
"""
@inline @fastmath function haversine(lon1, lat1, lon2, lat2)
    R = 6372.8  # Earth radius in kilometers

    dLat = deg2rad(lat2 - lat1)
    dLon = deg2rad(lon2 - lon1)
    lat1 = deg2rad(lat1)
    lat2 = deg2rad(lat2)
    lon1 = deg2rad(lon1)
    lon2 = deg2rad(lon2)
    a = sin(dLat/2)^2 + cos(lat1)*cos(lat2)*sin(dLon/2)^2
    c = 2*asin(sqrt(a))
    theta = atan(sin(dLon)*cos(lat2),
                 cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(dLon))
    theta = (rad2deg(theta) + 360) % 360
    return R*c*0.5399565, theta
end


function rotate_point(ox, oy, px, py, angle,ret_x)
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    qx = ox + cos(angle) * (px - ox) - sin(angle) * (py - oy)
    qy = oy + sin(angle) * (px - ox) + cos(angle) * (py - oy)
    if ret_x == true
        return qx
    else
        return qy
    end
end


function generate_grid(start_lon, start_lat, finish_lon, finish_lat, nodes)
    dist = haversine(start_lon, start_lat, finish_lon, finish_lat)
    spacing = dist[1]/(nodes+1)
    alpha = dist[2]
    x_dist = spacing
    y_dist = spacing
    grid_x = reshape(start_lon.+[i*x_dist for i in range(1, length=nodes) for j in range(0, length=nodes).-(nodes-1)/2], (nodes, nodes))
    grid_y = reshape(start_lat.+[j*y_dist for i in range(1, length=nodes) for j in range(0, length=nodes).-(nodes-1)/2], (nodes, nodes))
    rot_grid_x = [rotate_point(start_lon, start_lat, x, y, alpha, true) for (x, y) in zip(grid_x,grid_y)]
    rot_grid_y = [rotate_point(start_lon, start_lat, x, y, alpha, false) for (x, y) in zip(grid_x,grid_y)]
    return rot_grid_x, rot_grid_y
end


"""
    euclidean(x1, y1, x2, y2)

Calculate the distance and argument between the vector and north between points 1 and 2.
euclidean(0.0, 0.0, 10.0, 10.0) = (14.142135623730951, 45.0)
euclidean(0.0, 0.0, 10.0, 10.0) = (14.142135623730951, 45.0)
euclidean(0.0, 0.0, 10.0, 10.0) = (14.142135623730951, 45.0)
euclidean(0.0, 0.0, 10.0, 10.0) = (14.142135623730951, 45.0)
"""
@inline @fastmath function euclidean(x1::Float64, y1::Float64, x2::Float64, y2::Float64)
    dx = x2 - x1
    dy = y2 - y1
    dist = (dx^2 + dy^2)^(0.5)
    theta = 0.0
    if dy > 0.0 && dx > 0.0
        theta = 90.0 - rad2deg(atan(dy, dx))
    elseif dy == 0 && dx > 0.0
        theta = 90.0
    elseif dx == 0 && dy < 0.0
        theta = 180.0
    elseif dy < 0.0 && dx > 0.0
        theta = 90.0 + rad2deg(atan(abs(dy), dx))
    elseif dy < 0.0 && dx < 0.0
        theta = 180 + rad2deg(atan(abs(dy), abs(dx)))
    elseif dy > 0.0 && dx < 0.0
        theta = 360 - rad2deg(atan(dy, abs(dx)))
    end
    return dist, theta
end


"""Calculate the number of nodes for a specific distance in nm."""
function calc_nodes(lon1, lon2, lat1, lat2, req_d)
    d = haversine(lon1, lon2, lat1, lat2)[1]
    return round(Int, d/req_d)
end


