# Piecewise-linear reconstruction: Gouraud-shaded Delaunay triangles.
#
# Why this rather than the Voronoi cells of voronoibrot.jl. Nearest-neighbour
# reconstruction is piecewise *constant*, so its error goes like h·|∇ν| in the
# cell size h — halving the visible error costs four times the samples, and it
# costs them everywhere in the smooth region. A piecewise-*linear*
# reconstruction has error h²·|∇²ν|, so a region where the field is nearly
# linear needs almost no samples at all, while near ∂M every derivative blows
# up and you still need many. In other words the Voronoi is the reconstruction
# that makes *least* use of smoothness, and moving to triangles makes the
# "concentrate samples on the boundary" property stronger rather than weaker.
# Measured: this draws a smoother exterior from 2 300 samples than the Voronoi
# manages from 15 700.
#
# Two things are load-bearing:
#
#  * [`smoothcount`](@ref) is a prerequisite, not a nicety. The raw integer
#    count is a staircase, and interpolating a staircase linearly just gives
#    tilted stairs — facets aligned to the escape-time bands.
#  * Linear interpolation lies across ∂M. A triangle with both an interior and
#    an escaped vertex straddles the discontinuity, and Gouraud shading it
#    would smear the boundary across a whole cell. Those triangles are instead
#    subdivided into their three barycentric quads and flat-shaded, one per
#    vertex — which is exactly nearest-neighbour reconstruction, locally. So
#    the picture is linear where the field is smooth and nearest-neighbour
#    where it jumps, and the boundary lands halfway between samples instead of
#    being dilated to the nearest triangle edge.
#
# That composes with the refinement in refine.jl, which is already driving
# samples onto exactly those straddling edges: the flat-shaded band shrinks
# toward sub-pixel as you refine.
#
# The honest cost, for the record: a flat Voronoi cell says "everything here
# is nearest to a sample worth ν and I know nothing finer", which is true. A
# Gouraud triangle *asserts* that the field is linear across it, which is only
# approximately true. The half-finished picture stops being a strict lower
# bound.

"""
    MeshBrot

Handle returned by [`meshbrotplot!`](@ref). Same contract as
[`VoronoiBrot`](@ref) — the caller owns the `SampleSet` and calls
[`commit!`](@ref) after advancing it.
"""
mutable struct MeshBrot
    set::SampleSet
    axis::Any
    verts::Observable{Vector{Point2f}}
    faces::Observable{Vector{TriFace}}
    vcols::Observable{Vector{Float32}}
    colorrange::Observable{Vec2f}
    admit::Symbol
    smooth::Bool
    positions::Vector{Point2f}     # last triangulated generator set
    tri::Any                       # its Delaunay triangulation
    outside_color::Any
    attributes::Dict{Symbol, Any}
    plot::Any
    backdrop::Any
end

"""
    meshbrotplot!(ax, set::SampleSet; kwargs...) -> MeshBrot

Draw `set` as a Gouraud-shaded Delaunay mesh. Keywords match
[`voronoibrotplot!`](@ref) except that there is no `clip`: a mesh only covers
the convex hull of its samples, so there is nothing unbounded to clip. That
is also its one visible artefact — outside the hull nothing is drawn, so the
corners of a rectangular viewport show `outside_color` until samples reach
them. `sprinkle!` over the viewport fixes it; a bare cluster of clicks will
not.
"""
function meshbrotplot!(ax, set::SampleSet;
        admit = :immediate,
        smooth = true,
        colormap = :PRGn_9,
        colorscale = log10,
        interior_color = :black,
        pending_color = :grey30,
        outside_color = :black,
        kwargs...)

    admit in (:onresolve, :immediate) ||
        throw(ArgumentError("admit must be :onresolve or :immediate, got $admit"))

    verts = Observable(Point2f[])
    faces = Observable(TriFace[])
    vcols = Observable(Float32[])
    colorrange = Observable(Vec2f(1, 2))

    attributes = Dict{Symbol, Any}(
        :color => vcols,
        :colorrange => colorrange,
        :colormap => colormap,
        :colorscale => colorscale,
        :nan_color => interior_color,
        :lowclip => pending_color,
        :shading => NoShading,
    )
    merge!(attributes, kwargs)

    mb = MeshBrot(set, ax, verts, faces, vcols, colorrange, admit, smooth,
                  Point2f[], nothing, outside_color, attributes, nothing, nothing)
    commit!(mb)
    return mb
end

"""
    meshbrotplot(set::SampleSet; figure = (;), axis = (;), kwargs...)

Convenience wrapper returning `(fig, ax, mb)`.
"""
function meshbrotplot(set::SampleSet; figure = (;), axis = (;), kwargs...)
    fig = Figure(; figure...)
    ax = Axis(fig[1, 1]; aspect = DataAspect(), axis...)
    mb = meshbrotplot!(ax, set; kwargs...)
    return fig, ax, mb
end

# A vertex can be Gouraud-interpolated only if it carries a real value: not
# interior (no value at all) and not pending (no value *yet*).
_interpolatable(v) = isfinite(v) && v > PENDING_VALUE

"""
    build_mesh(positions, values, tri) -> (verts, faces, vcols)

Triangles whose three vertices all carry a real value are emitted as-is and
get Gouraud shading. The rest straddle a discontinuity — ∂M, or the edge of
what has been computed — and are split into three flat quads about the
centroid, one per vertex, each carrying that vertex's own colour. Every
midpoint and centroid is duplicated per owning vertex, precisely because the
three quads must *not* share a colour.
"""
function build_mesh(positions, values, tri)
    verts = copy(positions)
    vcols = copy(values)
    faces = TriFace[]
    for (a, b, c) in DelTri.each_solid_triangle(tri)
        if _interpolatable(values[a]) && _interpolatable(values[b]) &&
                _interpolatable(values[c])
            push!(faces, TriFace(a, b, c))
            continue
        end
        pa, pb, pc = positions[a], positions[b], positions[c]
        centroid = (pa + pb + pc) / 3
        mab, mbc, mca = (pa + pb) / 2, (pb + pc) / 2, (pc + pa) / 2
        for (v, m1, m2) in ((a, mab, mca), (b, mbc, mab), (c, mca, mbc))
            i = length(verts) + 1
            push!(verts, m1, centroid, m2)
            push!(vcols, values[v], values[v], values[v])
            push!(faces, TriFace(v, i, i + 1),
                         TriFace(v, i + 1, i + 2))
        end
    end
    return verts, faces, vcols
end

"""
    commit!(mb::MeshBrot)

Republish the sample set. The Delaunay triangulation is only rebuilt when the
generator positions actually change — resolving a sample changes its colour
and its triangle's *classification*, but not the triangulation, and swapping
a face list is a buffer update rather than an `O(n log n)` retriangulation.
"""
function commit!(mb::MeshBrot)
    positions, values, colorrange = sample_arrays(mb.set; admit = mb.admit,
                                                  smooth = mb.smooth)
    mb.colorrange[] = colorrange
    triangulable(positions) || return mb

    if positions != mb.positions
        mb.positions = positions
        mb.tri = DelTri.triangulate([(Float64(p[1]), Float64(p[2])) for p in positions];
                                    randomise = false)
    end

    verts, faces, vcols = build_mesh(positions, values, mb.tri)
    # colours first, then geometry: the lengths must never be observed apart
    mb.vcols[] = vcols
    mb.verts[] = verts
    mb.faces[] = faces

    if mb.plot === nothing
        # Backdrop first, so the region outside the convex hull reads as
        # "no data" rather than as whatever is behind the axis.
        mb.backdrop = poly!(mb.axis, lift(identity, mb.axis.finallimits);
                            color = mb.outside_color)
        translate!(mb.backdrop, 0, 0, -1)
        mb.plot = mesh!(mb.axis, mb.verts, mb.faces; mb.attributes...)
    end
    return mb
end

"""
    detach!(mb)

Remove the plot from its axis but keep the samples. `commit!` rebuilds it,
so this is how a caller swaps one reconstruction for another.
"""
function detach!(mb::MeshBrot)
    mb.plot === nothing || delete!(mb.axis, mb.plot)
    mb.backdrop === nothing || delete!(mb.axis, mb.backdrop)
    mb.plot = nothing
    mb.backdrop = nothing
    return mb
end

function Base.empty!(mb::MeshBrot)
    detach!(mb)
    empty!(mb.set)
    empty!(mb.positions)
    mb.verts[] = Point2f[]
    mb.faces[] = TriFace[]
    mb.vcols[] = Float32[]
    return mb
end
