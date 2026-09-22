# Shared plumbing for the reconstructions in voronoibrot.jl and meshbrot.jl.
#
# Both turn the same `SampleSet` into the same three things — generator
# positions, per-generator values, and a colour range — and differ only in
# what they draw with them. Voronoi is piecewise constant, the mesh is
# piecewise linear; see meshbrot.jl for why that matters more than it sounds.

# DelaunayTriangulation, reached through Makie's alias. It is a hard
# dependency of Makie (not a weakdep), so it is already in the manifest and
# BrotViz stays free of a direct dependency on it. Tripwire: if Makie ever
# drops the alias, add DelaunayTriangulation to Project.toml — it costs
# nothing to install.
const DelTri = Makie.DelTri
const TriFace = Makie.GeometryBasics.TriangleFace{Int}

# Sentinel value for a cell whose sample has not resolved yet. Below any real
# count (which is >= 1) so it routes to `lowclip`, and strictly positive so it
# survives the default `colorscale = log10`.
const PENDING_VALUE = 1.0f-6

"""
    triangulable(ps) -> Bool

Whether `ps` can be Delaunay-triangulated: at least three points, not all on
one line. The collinearity test uses the lexicographic extremes as a
baseline, so it is a single linear pass and scale-invariant.
"""
function triangulable(ps)
    length(ps) < 3 && return false
    a = argmin(p -> (p[1], p[2]), ps)
    b = argmax(p -> (p[1], p[2]), ps)
    dx, dy = Float64(b[1] - a[1]), Float64(b[2] - a[2])
    L2 = dx * dx + dy * dy
    L2 > 0 || return false
    for p in ps
        ex, ey = Float64(p[1] - a[1]), Float64(p[2] - a[2])
        # |cross| = |a-b| * distance(p, line ab); want that distance > 1e-6 |a-b|
        abs(dx * ey - dy * ex) > 1.0e-6 * L2 && return true
    end
    return false
end

"""
    limits_of(ax) -> ((xmin, xmax), (ymin, ymax))

The axis's current data-space viewport, as the samplers want it.
"""
function limits_of(ax)
    r = ax.finallimits[]
    o, w = r.origin, r.widths
    return ((o[1], o[1] + w[1]), (o[2], o[2] + w[2]))
end

"""
    sample_arrays(set; admit, smooth) -> (positions, values, colorrange)

One linear pass turning a `SampleSet` into what a reconstruction needs.

Building both arrays together is what guarantees the positional coupling
between a generator and its colour. `colorrange` spans the *escaped* values
only, so that interior (`NaN32`, routed to `nan_color`) and pending
(`PENDING_VALUE`, routed to `lowclip`) sit deliberately outside it. It never
starts below 1, because the default colorscale is `log10` and `log10(0)`
would poison the colour normalisation.
"""
function sample_arrays(set::SampleSet; admit::Symbol, smooth::Bool)
    immediate = admit === :immediate
    positions = Point2f[]
    values = Float32[]
    lo, hi = Inf32, -Inf32
    for s in set.samples
        if s.status === Pending
            immediate || continue
            push!(positions, Point2f(real(s.c), imag(s.c)))
            push!(values, PENDING_VALUE)
        else
            push!(positions, Point2f(real(s.c), imag(s.c)))
            v = Float32(value(set, s; smooth))
            push!(values, v)
            if isfinite(v)
                lo = min(lo, v)
                hi = max(hi, v)
            end
        end
    end
    colorrange = isfinite(lo) ? Vec2f(lo, hi > lo ? hi : lo + 1) : Vec2f(1, 2)
    return positions, values, colorrange
end
