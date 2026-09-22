# Annular wedge diagrams for digit expansions and kneading sequences.
#
# Both diagrams are built the same way: ring `n` partitions the circle of
# angles into the sectors on which the `n`-th symbol of θ is constant, and is
# drawn between radius `n-1` and `n`. Reading outward along the ray at angle θ
# therefore spells out the first `depth` symbols of θ.
#
# Base-`d` expansion. The `n`-th digit of θ is
#
#     x_n = floor(d^n θ) mod d,
#
# so ring `n` has `d^n` equal sectors with boundaries at multiples of `1/d^n`,
# and sector `k` carries digit `k mod d`. For d = 2 the innermost ring is a
# disk split into the upper half (digit 0) and lower half (digit 1), and each
# ring out from there doubles the sector count, alternating.
#
# Kneading sequence. `K(θ)_n` is the label of the arc of the critical
# partition containing `d^(n-1) θ`. That partition is cut by the `d` preimages
# `(θ+b)/d` of θ, so
#
#     d^(n-1) θ ∈ ((θ+b)/d, (θ+b+1)/d)  ⟺  (d^n - 1) θ mod d ∈ (b, b+1),
#
# giving
#
#     K(θ)_n = floor((d^n - 1) θ) mod d.
#
# Ring `n` therefore has `d^n - 1` equal sectors with boundaries at multiples
# of `1/(d^n - 1)`, and sector `k` carries symbol `k mod d`. The sector
# boundaries are precisely the angles whose `n`-th symbol is a star: they are
# the θ for which `(d^n - 1) θ` is an integer, i.e. for d = 2 the angles
# `k/(2^n - 1)` whose orbit lands on a ray through the critical point at step
# `n`. Stars are a measure-zero set, so they are drawn as lines rather than
# sectors (see `starwidth`).
#
# The two diagrams are the same construction with `d^n` replaced by `d^n - 1`,
# which is the whole difference between the fixed partition at {0, 1/2, ...}
# and the θ-dependent critical partition.

# Cap on the total number of wedge polygons, to catch a `depth` that would
# hang the backend rather than produce a figure.
const MAX_WEDGES = 100_000

function _wedgecolors(degree::Integer, colors)
    if colors === nothing
        return _defaultwedgecolors(degree)
    end
    cs = [Makie.to_color(c) for c in colors]
    length(cs) == degree ||
        throw(ArgumentError("degree $degree needs $degree colors, got $(length(cs))"))
    return cs
end

# Degree 2 gets the A/B colors used by the tree plots. Higher degrees fall back
# to Wong's categorical palette, then to a large categorical scheme.
function _defaultwedgecolors(degree::Integer)
    degree == 2 && return [Makie.to_color(:royalblue), Makie.to_color(:orangered)]
    wong = Makie.wong_colors()
    degree <= length(wong) && return collect(wong[1:degree])
    return [get(ColorSchemes.glasbey_hv_n256, (k - 1) / (degree - 1)) for k in 1:degree]
end

# Validated up front by both the mutating and non-mutating methods, so that a
# bad `depth` errors before an empty axis is built.
function _checkwedgeargs(depth::Integer, degree::Integer, nwedges::Integer)
    depth >= 1 || throw(ArgumentError("depth must be at least 1, got $depth"))
    degree >= 2 || throw(ArgumentError("degree must be at least 2, got $degree"))
    nwedges <= MAX_WEDGES || throw(ArgumentError(
        "depth=$depth at degree=$degree needs $nwedges wedges, over the limit of $MAX_WEDGES"))
    return nothing
end

# Counted in BigInt: degree^depth overflows Int64 well before it stops being a
# plausible thing to type.
_expansionwedgecount(depth::Integer, degree::Integer) =
    depth >= 1 && degree >= 2 ? sum(big(degree)^n for n in 1:depth) : big(0)

_kneadingwedgecount(depth::Integer, degree::Integer) =
    depth >= 1 && degree >= 2 ? sum(big(degree)^n - 1 for n in 1:depth) : big(0)

"""
    _ringwedges!(ax, nsectors, r0, r1, colors; strokecolor, strokewidth, arcstep)

Draw one ring of `nsectors` equal wedges between radii `r0` and `r1`, with
sector `k` (counted counterclockwise from θ = 0) colored `colors[k mod end + 1]`.
Expects `ax` to be a `PolarAxis`, so vertices are `Point2f(θ, r)`.
"""
function _ringwedges!(ax, nsectors::Integer, r0::Real, r1::Real, colors;
        strokecolor=:transparent, strokewidth=0.0, arcstep=0.05)
    ncolors = length(colors)
    span = 2π / nsectors
    # A whole ring of one sector still needs enough vertices to look round.
    nseg = max(2, ceil(Int, span / arcstep) + 1)
    for k in 0:(nsectors - 1)
        θs = range(span * k, span * (k + 1); length=nseg)
        wedge = Vector{Point2f}(undef, 2nseg)
        @inbounds for i in 1:nseg
            wedge[i] = Point2f(θs[i], r1)
            wedge[2nseg - i + 1] = Point2f(θs[i], r0)
        end
        poly!(ax, wedge; color=colors[mod(k, ncolors) + 1],
            strokecolor=strokecolor, strokewidth=strokewidth)
    end
    return ax
end

"""
    _ringstars!(ax, nsectors, r0, r1; color, width)

Draw the `nsectors` sector boundaries of a ring as radial lines. For the
kneading diagram these are exactly the angles whose symbol at this depth is a
star.
"""
function _ringstars!(ax, nsectors::Integer, r0::Real, r1::Real; color=:black, width=0.0)
    width > 0 || return ax
    span = 2π / nsectors
    for k in 0:(nsectors - 1)
        lines!(ax, [Point2f(span * k, r0), Point2f(span * k, r1)];
            color=color, linewidth=width)
    end
    return ax
end

function _wedgeaxis(fig, position, rmax)
    ax = PolarAxis(position;
        rticklabelsvisible=false,
        thetaticklabelsvisible=false,
        rgridvisible=false,
        thetagridvisible=false,
        spinevisible=false)
    hidedecorations!(ax)
    rlims!(ax, 0, rmax)
    return ax
end

# ============================================================================
# Digit expansions
# ============================================================================

"""
    expansionwedgeplot!(ax, depth; degree=2, colors=nothing, innerradius=0.0,
                        ringwidth=1.0, strokecolor=:gray60, strokewidth=0.5,
                        arcstep=0.05)

Draw the first `depth` rings of the base-`degree` digit expansion diagram into
an existing `PolarAxis`.

Ring `n` has `degree^n` equal sectors with boundaries at multiples of
`1/degree^n`; the sector starting at `k/degree^n` carries digit `k mod degree`,
because the `n`-th digit of θ is `floor(degree^n * θ) mod degree`. Travelling
outward along the ray at angle θ therefore reads off the first `depth` digits
of θ. For `degree = 2` the innermost ring is a disk whose upper half is the
digit-0 color and whose lower half is the digit-1 color.

`colors` must have one entry per digit, in the order `0, 1, ..., degree-1`; it
defaults to a palette chosen for the degree. Rings are drawn between radii
`innerradius + (n-1)*ringwidth` and `innerradius + n*ringwidth`. `arcstep` is
the maximum angular step, in radians, used to polygonize each arc.

See also [`kneadingwedgeplot!`](@ref).
"""
function expansionwedgeplot!(ax, depth::Integer; degree::Integer=2, colors=nothing,
        innerradius::Real=0.0, ringwidth::Real=1.0,
        strokecolor=:gray60, strokewidth=0.5, arcstep::Real=0.05)
    _checkwedgeargs(depth, degree, _expansionwedgecount(depth, degree))
    cs = _wedgecolors(degree, colors)

    for n in 1:depth
        _ringwedges!(ax, degree^n,
            innerradius + (n - 1) * ringwidth, innerradius + n * ringwidth, cs;
            strokecolor=strokecolor, strokewidth=strokewidth, arcstep=arcstep)
    end
    return ax
end

"""
    expansionwedgeplot(depth; size=(600, 600), kwargs...)

Digit expansion wedge diagram in a new `Figure` with a bare `PolarAxis`.
Returns `(fig, ax)`. See [`expansionwedgeplot!`](@ref) for the keyword
arguments.

```julia
fig, ax = expansionwedgeplot(6)
fig, ax = expansionwedgeplot(4; colors=(:white, :black))
fig, ax = expansionwedgeplot(4; degree=3)
```
"""
function expansionwedgeplot(depth::Integer; degree::Integer=2, size=(600, 600),
        innerradius::Real=0.0, ringwidth::Real=1.0, kwargs...)
    _checkwedgeargs(depth, degree, _expansionwedgecount(depth, degree))
    fig = Figure(size=size)
    ax = _wedgeaxis(fig, fig[1, 1], innerradius + depth * ringwidth)
    expansionwedgeplot!(ax, depth; degree=degree,
        innerradius=innerradius, ringwidth=ringwidth, kwargs...)
    return fig, ax
end

# ============================================================================
# Kneading sequences
# ============================================================================

"""
    kneadingwedgeplot!(ax, depth; degree=2, colors=nothing, innerradius=0.0,
                       ringwidth=1.0, starcolor=:black, starwidth=0.0,
                       strokecolor=:transparent, strokewidth=0.0, arcstep=0.05)

Draw the first `depth` rings of the kneading sequence diagram into an existing
`PolarAxis`.

Ring `n` has `degree^n - 1` equal sectors with boundaries at multiples of
`1/(degree^n - 1)`; the sector starting at `k/(degree^n - 1)` carries symbol
`k mod degree`, because

    K(θ)_n = floor((degree^n - 1) * θ) mod degree.

Travelling outward along the ray at angle θ therefore reads off the first
`depth` symbols of the kneading sequence of θ, with symbol `0` being `A`. This
is the same construction as [`expansionwedgeplot!`](@ref) with `degree^n`
replaced by `degree^n - 1`, which is the entire difference between the fixed
partition of the circle and the θ-dependent critical partition.

The sector boundaries are exactly the angles whose `n`-th symbol is a star, a
measure-zero set. Set `starwidth` to a positive number to draw them as radial
lines in `starcolor`; this is legible at small `depth` and turns the outer
rings into a solid band at large `depth`, so it is off by default.

`colors` must have one entry per symbol, in the order `A, B, ...`; it defaults
to a palette chosen for the degree.
"""
function kneadingwedgeplot!(ax, depth::Integer; degree::Integer=2, colors=nothing,
        innerradius::Real=0.0, ringwidth::Real=1.0,
        starcolor=:black, starwidth::Real=0.0,
        strokecolor=:transparent, strokewidth=0.0, arcstep::Real=0.05)
    _checkwedgeargs(depth, degree, _kneadingwedgecount(depth, degree))
    cs = _wedgecolors(degree, colors)

    for n in 1:depth
        r0 = innerradius + (n - 1) * ringwidth
        r1 = innerradius + n * ringwidth
        nsectors = degree^n - 1
        _ringwedges!(ax, nsectors, r0, r1, cs;
            strokecolor=strokecolor, strokewidth=strokewidth, arcstep=arcstep)
        _ringstars!(ax, nsectors, r0, r1; color=starcolor, width=starwidth)
    end
    return ax
end

"""
    kneadingwedgeplot(depth; size=(600, 600), kwargs...)

Kneading sequence wedge diagram in a new `Figure` with a bare `PolarAxis`.
Returns `(fig, ax)`. See [`kneadingwedgeplot!`](@ref) for the keyword
arguments.

```julia
fig, ax = kneadingwedgeplot(6)
fig, ax = kneadingwedgeplot(4; starwidth=1.0)
fig, ax = kneadingwedgeplot(4; degree=3)
```
"""
function kneadingwedgeplot(depth::Integer; degree::Integer=2, size=(600, 600),
        innerradius::Real=0.0, ringwidth::Real=1.0, kwargs...)
    _checkwedgeargs(depth, degree, _kneadingwedgecount(depth, degree))
    fig = Figure(size=size)
    ax = _wedgeaxis(fig, fig[1, 1], innerradius + depth * ringwidth)
    kneadingwedgeplot!(ax, depth; degree=degree,
        innerradius=innerradius, ringwidth=ringwidth, kwargs...)
    return fig, ax
end

# ============================================================================
# Symbol lookup
# ============================================================================

"""
    expansiondigits(theta, depth; degree=2)

The first `depth` base-`degree` digits of the angle `theta`, as read off the
diagram drawn by [`expansionwedgeplot!`](@ref).
"""
function expansiondigits(theta, depth::Integer; degree::Integer=2)
    t = mod(theta, 1)
    return [Int(mod(floor(BigInt, big(degree)^n * t), degree)) for n in 1:depth]
end

"""
    kneadingsymbols(theta, depth; degree=2)

The first `depth` kneading symbols of the angle `theta` as symbol indices
`0, 1, ..., degree-1` (`0` being `A`), as read off the diagram drawn by
[`kneadingwedgeplot!`](@ref). Returns `nothing` in place of an index where
`theta` lands exactly on a sector boundary, i.e. where the symbol is a star.
"""
function kneadingsymbols(theta, depth::Integer; degree::Integer=2)
    t = mod(theta, 1)
    out = Union{Int, Nothing}[]
    for n in 1:depth
        x = (big(degree)^n - 1) * t
        push!(out, isinteger(x) ? nothing : Int(mod(floor(BigInt, x), degree)))
    end
    return out
end
