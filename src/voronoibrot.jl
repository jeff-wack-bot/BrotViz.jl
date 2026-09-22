# Voronoi reconstruction of a scattered SampleSet.
#
# Makie-only — no GLMakie — so this renders under CairoMakie too and can
# produce static figures. Follows the plain-function style of recipes.jl
# (`foo!(ax, …)` plus `foo(…)`), not `@recipe`.
#
# Each resolved sample paints the cell of points nearest to it: a flat cell
# says "everything here is nearest to a sample worth ν, and I know nothing
# finer". Because the generators are complex numbers and the cells are
# polygons in data space, pan and zoom are exact at every scale and old data
# never goes stale.
#
# Gotchas verified against Makie 0.24.15 / DelaunayTriangulation 1.6.7:
#
#  * `triangulate` throws on fewer than three points *and* on collinear point
#    sets, and Makie's compute pipeline resolves the tessellation eagerly as
#    soon as either the positions or the colors are written. `visible = false`
#    does **not** protect you — the plot must never exist while its generators
#    are degenerate. So `voronoiplot!` is called lazily, on the first commit
#    whose generator set is triangulable, and the plot is deleted again if the
#    set is emptied.
#  * `voronoiplot` asserts `length(color) == num_points(triangulation)` and
#    indexes `color` by *generator index*, so `gens` and `vals` are
#    length-coupled and positionally coupled. `commit!` rebuilds both from
#    scratch in the same pass, which makes that true by construction.

"""
    VoronoiBrot

Handle returned by [`voronoibrotplot!`](@ref). Owns the observables that
back the plot; the caller owns the `SampleSet` and calls [`commit!`](@ref)
whenever it has advanced it.

`admit` selects the cell admission policy and can be flipped at any time:

  * `:immediate` (default) — every sample gets a cell as soon as it is
    queued, pending ones in `pending_color`. A batch of new samples tiles
    the plane grey at once and then colours in one cell at a time as the
    scheduler resolves them.
  * `:onresolve` — only resolved samples get a cell, so a pending sample is
    simply not drawn.

`:immediate` is the default for a reason beyond looks. `voronoiplot`
retriangulates whenever the *generator positions* change but not when only
the colours do, and at 16 000 samples that is 240 ms against 15 ms. Under
`:onresolve` every single sample that resolves changes the generator set, so
the whole tessellation is rebuilt on every frame that makes progress —
about 4 fps. Under `:immediate` the generators are fixed for the lifetime of
a batch: one retessellation when the batch is queued, then colour-only
updates as it resolves. Same data, same recipe, 16× cheaper.
"""
mutable struct VoronoiBrot
    set::SampleSet
    axis::Any
    gens::Observable{Vector{Point2f}}
    vals::Observable{Vector{Float32}}
    colorrange::Observable{Vec2f}
    admit::Symbol
    smooth::Bool
    attributes::Dict{Symbol, Any}
    plot::Any                      # the voronoiplot, or nothing until triangulable
end

"""
    voronoibrotplot!(ax, set::SampleSet; kwargs...) -> VoronoiBrot

Draw `set` into `ax` as a Voronoi tessellation. Returns a handle to pass to
`commit!`.

Keywords:

  * `admit = :immediate` — cell admission policy, see [`VoronoiBrot`](@ref).
  * `smooth = true` — colour by the renormalised count ([`smoothcount`](@ref))
    rather than the raw integer, which removes the escape-time banding.
  * `colormap = :PRGn_9` — escape-count colormap, as the raster renderer uses.
  * `colorscale = log10` — escape counts are wildly skewed: over the classic
    view most exterior samples leave in under ten iterations while a few near
    ∂M take hundreds, so a linear ramp collapses the whole exterior onto one
    colour. A log scale spreads the equipotential bands out. Pass `identity`
    for the raw count.
  * `interior_color = :black` — non-escaping cells. Routed through `nan_color`,
    the same convention as `mandelbrotsetplot`.
  * `pending_color = :grey30` — cells of samples that have not resolved yet.
  * `clip = :limits` — clip cells to `ax.finallimits`, so the unbounded outer
    cells do not sprawl. Pass a `Rect2`/`Circle` for a fixed clip, or
    `Makie.automatic` for Makie's own bounds. Anything other than `:limits`
    is handed to `voronoiplot` unchanged.

Any further keywords are forwarded to `voronoiplot!`.
"""
function voronoibrotplot!(ax, set::SampleSet;
        admit = :immediate,
        smooth = true,
        colormap = :PRGn_9,
        colorscale = log10,
        interior_color = :black,
        pending_color = :grey30,
        strokewidth = 0,
        clip = :limits,
        kwargs...)

    admit in (:onresolve, :immediate) ||
        throw(ArgumentError("admit must be :onresolve or :immediate, got $admit"))

    gens = Observable(Point2f[])
    vals = Observable(Float32[])
    # Starts at (1, 2), not (0, 1): the default colorscale is log10 and
    # log10(0) = -Inf poisons the colour normalisation. Real escape counts
    # are >= 1, so this range is also what `commit!` will only ever widen.
    colorrange = Observable(Vec2f(1, 2))

    attributes = Dict{Symbol, Any}(
        :color => vals,
        :colorrange => colorrange,
        :colormap => colormap,
        :colorscale => colorscale,
        :nan_color => interior_color,
        :lowclip => pending_color,
        :strokewidth => strokewidth,
        :show_generators => false,
        :clip => clip === :limits ? lift(identity, ax.finallimits) : clip,
    )
    merge!(attributes, kwargs)

    vb = VoronoiBrot(set, ax, gens, vals, colorrange, admit, smooth, attributes, nothing)

    commit!(vb)
    return vb
end

"""
    voronoibrotplot(set::SampleSet; figure = (;), axis = (;), kwargs...)

Convenience wrapper returning `(fig, ax, vb)`.
"""
function voronoibrotplot(set::SampleSet; figure = (;), axis = (;), kwargs...)
    fig = Figure(; figure...)
    ax = Axis(fig[1, 1]; aspect = DataAspect(), axis...)
    vb = voronoibrotplot!(ax, set; kwargs...)
    return fig, ax, vb
end

"""
    commit!(vb::VoronoiBrot)

Republish the sample set to the plot.

Rebuilding both arrays in one pass is what guarantees the length- and
order-coupling that `voronoiplot` asserts on. The generators are only
*written* if they actually changed: an equal-but-freshly-allocated vector
still counts as a change to Makie's compute pipeline and would cost a full
retriangulation, which is the difference between 240 ms and 15 ms at 16 000
samples. Comparing the two vectors first is `O(n)` and costs microseconds.
"""
function commit!(vb::VoronoiBrot)
    gens, vals, colorrange = sample_arrays(vb.set; admit = vb.admit, smooth = vb.smooth)
    vb.colorrange[] = colorrange

    if triangulable(gens)
        vb.vals[] = vals
        # Only touch the generators when they really moved — see above.
        gens == vb.gens[] || (vb.gens[] = gens)
        if vb.plot === nothing
            vb.plot = voronoiplot!(vb.axis, vb.gens; vb.attributes...)
        end
    end
    return vb
end


"""
    detach!(vb)

Remove the plot from its axis but keep the samples. `commit!` rebuilds it,
so this is how a caller swaps one reconstruction for another.
"""
function detach!(vb::VoronoiBrot)
    vb.plot === nothing || delete!(vb.axis, vb.plot)
    vb.plot = nothing
    return vb
end

"""
    Base.empty!(vb::VoronoiBrot)

Drop every sample and tear the tessellation back down. The plot has to go
before the generators do, since `voronoiplot` cannot survive being emptied.
"""
function Base.empty!(vb::VoronoiBrot)
    detach!(vb)
    empty!(vb.set)
    vb.gens[] = Point2f[]
    vb.vals[] = Float32[]
    return vb
end
