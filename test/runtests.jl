using MeshToSDF
using GeometryBasics
using Meshing
using Test

function sd_sphere(point, s)
    return √(sum(abs2, point)) - s
end

@testset "MeshToSDF.jl" begin
    @testset "Basic sphere test" begin
       # generate the sdf of a sphere
        n = 128
        level = 2 / n
        mesh_scale = 0.8
        r = range(-1, 1, n)
        sphere_sdf = [sd_sphere(Point3(x, y, z), mesh_scale) for x in r, y in r, z in r]

        # now generate a mesh out of it using marching cubes
        sdfv = sphere_sdf
        _vertices, _faces = isosurface(sdfv)
        _vertices = convert.(Point3{Float32}, _vertices)

        # now apply the level-set algorithm
        sdf_out = MeshToSDF.signed_distance_field(_vertices, _faces, n)

        # test of the maximum error is less than the expected resolution
        @test maximum(abs.(sdf_out) .- abs.(sphere_sdf)) < 2*level
    end
end
