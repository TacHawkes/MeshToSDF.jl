mag2(x) = sum(x->x^2, x)
function point_segment_distance(x0, x1, x2)
    dx = x2 - x1
    m2 = mag2(dx)
    # find parameter value of closest point on segment
    s12 = dot(x2 - x0, dx) / m2
    if s12 < 0
        s12 = zero(s12)
    elseif s12 > 1
        s12 = one(s12)
    end

    return norm(x0 - (s12*x1 + (1 - s12)*x2))
end

function point_triangle_distance(x0::P, x1::P, x2::P, x3::P) where {T <: Real, P <: Point3{T}}
    x13 = x1 - x3
    x23 = x2 - x3
    x03 = x0 - x3

    m13, m23, d = mag2(x13), mag2(x23), dot(x13, x23)
    invdet = one(T) / max(m13*m23 - d*d, T(1e-30))
    a, b = dot(x13, x03), dot(x23, x03)

    w23 = invdet * (m23*a - d*b)
    w31 = invdet * (m13*b - d*a)
    w12 = one(T) - w23 - w31
    if (w23 >= zero(T) && w31 >= zero(T) && w12 >= zero(T))
        # inside point_triangle_distance
        return norm(x0 - (w23*x1 + w31*x2 + w12*x3))
    else
        if (w23 > zero(T))
            return min(point_segment_distance(x0, x1, x2), point_segment_distance(x0, x1, x3))
        elseif (w31 > zero(T))
            return min(point_segment_distance(x0, x1, x2), point_segment_distance(x0, x2, x3))
        else
            return min(point_segment_distance(x0, x1, x3), point_segment_distance(x0, x2, x3))
        end
    end
end

function check_neighbour(tri, x, phi, closest_tri, gx, i0, j0, k0, i1, j1, k1)
    if closest_tri[i1, j1, k1] >= 1
        p, q, r = tri[closest_tri[i1, j1, k1]]
        d = point_triangle_distance(gx, x[p], x[q], x[r])
        if d < phi[i0, j0, k0]
            phi[i0, j0, k0] = d
            closest_tri[i0, j0, k0] = closest_tri[i1, j1, k1]
        end
    end

    return nothing
end

function sweep(tri, x, phi, closest_tri, origin, dx, di, dj, dk)
    # Determine the sweep order based on di, dj, dk
    i0, i1 = di > 0 ? (2, size(phi, 1)) : (size(phi, 1) - 1, 1)
    j0, j1 = dj > 0 ? (2, size(phi, 2)) : (size(phi, 2) - 1, 1)
    k0, k1 = dk > 0 ? (2, size(phi, 3)) : (size(phi, 3) - 1, 1)

    # Perform the sweep
    for k in (dk > 0 ? (k0:k1) : (k0:-1:k1)), j in (dj > 0 ? (j0:j1) : (j0:-1:j1)), i in (di > 0 ? (i0:i1) : (i0:-1:i1))
        gx = Point3((i-1) * dx + origin[1], (j-1) * dx + origin[2], (k-1) * dx + origin[3])
        check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i - di, j, k)
        check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i, j - dj, k)
        check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i - di, j - dj, k)
        check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i, j, k - dk)
        check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i - di, j, k - dk)
        check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i, j - dj, k - dk)
        check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i - di, j - dj, k - dk)
    end
end

function orientation(x1, y1, x2, y2)
    twice_signed_area = y1*x2 - x1*y2

    sign = begin
        if twice_signed_area > 0
            1
        elseif twice_signed_area < 0
            -1
        elseif y2 > y1
            1
        elseif y2 < y1
            -1
        elseif x1 > x2
            1
        elseif x1 < x2
            -1
        else
            0
        end
    end

    return sign, twice_signed_area
end

point_in_triangle_2d(x0, y0, x1, y1, x2, y2, x3, y3) = point_in_triangle_2d(promote(x0, y0, x1, y1, x2, y2, x3, y3)...)
function point_in_triangle_2d(x0::T, y0::T, x1::T, y1::T, x2::T, y2::T, x3::T, y3::T) where T <: Real
    x1 -= x0; x2 -= x0; x3 -= x0;
    y1 -= y0; y2 -= y0; y3 -= y0;
    signa, a = orientation(x2, y2, x3, y3)
    signa == 0 && return false, zero(T), zero(T), zero(T)
    signb, b = orientation(x3, y3, x1, y1)
    signb != signa && return false, zero(T), zero(T), zero(T)
    signc, c = orientation(x1, y1, x2, y2)
    signc != signa && return false, zero(T), zero(T), zero(T)
    sum = a + b + c

    @assert sum != 0

    a /= sum
    b /= sum
    c /= sum

    return true, a, b, c
end

function make_level_set3(tri, x, origin, dx, ni, nj, nk, exact_band=1)
    phi = fill((ni + nj + nk)*dx, (ni, nj, nk))
    closest_tri = zeros(Int, ni, nj, nk)
    intersection_count = zeros(Int, ni, nj, nk)
    # Initialize distances near the mesh and figure out intersection counts
    for t in eachindex(tri)
        p, q, r = tri[t]
        fip, fjp, fkp = ((convert(Point3{Float64}, x[p]) - origin) / dx) .+ 1
        fiq, fjq, fkq = ((convert(Point3{Float64}, x[q]) - origin) / dx) .+ 1
        fir, fjr, fkr = ((convert(Point3{Float64}, x[r]) - origin) / dx) .+ 1

        i0 = clamp(trunc(Int, min(fip, fiq, fir)) - exact_band, 1, ni)
        i1 = clamp(trunc(Int, max(fip, fiq, fir)) + exact_band + 1, 1, ni)
        j0 = clamp(trunc(Int, min(fjp, fjq, fjr)) - exact_band, 1, nj)
        j1 = clamp(trunc(Int, max(fjp, fjq, fjr)) + exact_band + 1, 1, nj)
        k0 = clamp(trunc(Int, min(fkp, fkq, fkr)) - exact_band, 1, nk)
        k1 = clamp(trunc(Int, max(fkp, fkq, fkr)) + exact_band + 1, 1, nk)

        for k in k0:k1, j in j0:j1, i in i0:i1
            gx = Point3((i-1) * dx + origin[1], (j-1) * dx + origin[2], (k-1) * dx + origin[3])
            d = point_triangle_distance(gx, x[p], x[q], x[r])
            if d < phi[i, j, k]
                phi[i, j, k] = d
                closest_tri[i, j, k] = t
            end
        end

        # And do intersection counts
        j0 = clamp(ceil(Int, min(fjp, fjq, fjr)), 1, nj)
        j1 = clamp(floor(Int, max(fjp, fjq, fjr)), 1, nj)
        k0 = clamp(ceil(Int, min(fkp, fkq, fkr)), 1, nk)
        k1 = clamp(floor(Int, max(fkp, fkq, fkr)), 1, nk)
        for k in k0:k1, j in j0:j1
            pit, a, b, c = point_in_triangle_2d(j, k, fjp, fkp, fjq, fkq, fjr, fkr)
            if pit
                fi = a * fip + b * fiq + c * fir  # intersection i coordinate
                i_interval = ceil(Int, fi)  # intersection is in (i_interval-1, i_interval]
                if i_interval < 1
                    intersection_count[1, j, k] += 1  # enlarge the first interval
                elseif i_interval <= ni
                    intersection_count[i_interval, j, k] += 1
                end
                # Ignoring intersections beyond the +x side of the grid
            end
        end
    end

    # Fill in the rest of the distances with fast sweeping
    for pass in 1:2
        sweep(tri, x, phi, closest_tri, origin, dx, +1, +1, +1)
        sweep(tri, x, phi, closest_tri, origin, dx, -1, -1, -1)
        sweep(tri, x, phi, closest_tri, origin, dx, +1, +1, -1)
        sweep(tri, x, phi, closest_tri, origin, dx, -1, -1, +1)
        sweep(tri, x, phi, closest_tri, origin, dx, +1, -1, +1)
        sweep(tri, x, phi, closest_tri, origin, dx, -1, +1, -1)
        sweep(tri, x, phi, closest_tri, origin, dx, +1, -1, -1)
        sweep(tri, x, phi, closest_tri, origin, dx, -1, +1, +1)
    end

    # Determine signs (inside/outside) from intersection counts
    for k in 1:nk, j in 1:nj
        total_count = 0
        for i in 1:ni
            total_count += intersection_count[i, j, k]
            if isodd(total_count)
                phi[i, j, k] = -phi[i, j, k]  # inside the mesh
            end
        end
    end

    return phi
end

function signed_distance_field(V, F, size)
    _bbmin = Point3{Float32}(-1.0, -1.0, -1.0)
    dx = 2.0 / size
    return make_level_set3(F, V, _bbmin, Float32(dx), size, size, size)
end


function normalize_mesh(V, mesh_scale=0.8)
    bbmin = Point3(minimum(p[1] for p in V), minimum(p[2] for p in V), minimum(p[3] for p in V))
    bbmax = Point3(maximum(p[1] for p in V), maximum(p[2] for p in V), maximum(p[3] for p in V))
    center = (bbmin + bbmax) * 0.5
    scale = 2.0 * mesh_scale / maximum(bbmax - bbmin)
    return (V .- center) * scale
end
