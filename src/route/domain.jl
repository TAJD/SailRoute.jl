using PyCall

machinery = pyimport("importlib.machinery")
loader = machinery.SourceFileLoader("pydomain",ENV["HOME"]*"/sail_route_old/src/route/pydomain.py")
pd = loader.load_module("pydomain")


"""Calculate the minimum distance between two angles."""
@inline @fastmath function min_angle(a, b)
    abs(mod(a - b + 180.0, 360.0) - 180.0)
end


"""
    haversine(lon1, lat1, lon2, lat2)

Calculate the haversine distance and bearing. Distance is in nm.
"""
@inline @fastmath function haversine(lon1::Float64, lat1::Float64, lon2::Float64, lat2::Float64)
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
    req_n = round(Int, d/req_d)
    req_n
end


