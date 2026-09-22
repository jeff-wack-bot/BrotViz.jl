# Resumable escape-time samples for the scatter/Voronoi explorer.
#
# Pure compute: no Makie here. The contract is strict round-robin —
# `run_round!` advances *every* pending sample by the same number of
# iterations, so at any instant all live samples have received the same
# amount of work. That is what makes a half-finished picture a complete
# lower bound rather than a half-drawn one.
#
# This is deliberately a parallel implementation of `EscapeTimeProblem` /
# `escapetime` (renderfractal.jl) rather than a replacement. That API runs to
# completion and its `f`/`stop` fields are `::Function`, which costs a dynamic
# dispatch on every iteration; here the kernel is concrete and resumable.
# `test/runtests.jl` holds the regression test that keeps the two in step.

@enum SampleStatus Pending Escaped Interior

mutable struct Sample
    c::ComplexF64      # parameter — also the sample location in the plane
    z::ComplexF64      # current orbit position
    n::Int             # iterations performed so far
    status::SampleStatus
end

Sample(c) = Sample(ComplexF64(c), 0.0im, 0, Pending)

"""
    SampleSet(; escape_radius2 = 100.0, maxiter = 512)

A growable bag of samples plus the indices of those still iterating.

`escape_radius2` is the *squared* escape radius, matching the convention of
`escape(R)` in renderfractal.jl (which tests `abs2(z) > R`). The default of
`100.0` is the `escape(100)` used by the raster renderer, i.e. |z| > 10.
"""
struct SampleSet
    samples::Vector{Sample}
    pending::Vector{Int}    # indices into `samples`, compacted each round
    escape_radius2::Float64
    maxiter::Int
end

SampleSet(; escape_radius2 = 100.0, maxiter = 512) =
    SampleSet(Sample[], Int[], Float64(escape_radius2), Int(maxiter))

Base.length(set::SampleSet) = length(set.samples)
Base.isempty(set::SampleSet) = isempty(set.samples)

npending(set::SampleSet) = length(set.pending)
nresolved(set::SampleSet) = length(set.samples) - length(set.pending)

function Base.empty!(set::SampleSet)
    empty!(set.samples)
    empty!(set.pending)
    return set
end

"""
    push_sample!(set, c; mindist = 0.0)

Enqueue a new pending sample at parameter `c` and return its index.

If `mindist > 0` and `c` lies within that distance of an existing sample,
nothing is added and `nothing` is returned. Coincident generators make
DelaunayTriangulation warn and drop points, so callers that take positions
from mouse clicks should pass a tolerance scaled to the current viewport.
The check is a linear scan, which is fine at click cadence.
"""
function push_sample!(set::SampleSet, c; mindist = 0.0)
    cc = ComplexF64(c)
    if mindist > 0
        d2 = mindist * mindist
        for s in set.samples
            abs2(s.c - cc) < d2 && return nothing
        end
    end
    push!(set.samples, Sample(cc))
    push!(set.pending, length(set.samples))
    return length(set.samples)
end

"""
    advance!(s, chunk, R2, maxiter) -> Bool

Iterate `z ↦ z² + c` at most `chunk` more times, stopping early if `s`
resolves. Returns `true` if `s` is now resolved (`Escaped` or `Interior`).

The hot loop: concretely typed and non-allocating.
"""
function advance!(s::Sample, chunk, R2, maxiter)
    z, c, n = s.z, s.c, s.n
    for _ in 1:chunk
        z = z * z + c
        n += 1
        if abs2(z) > R2
            s.z, s.n, s.status = z, n, Escaped
            return true
        elseif n >= maxiter
            s.z, s.n, s.status = z, n, Interior
            return true
        end
    end
    s.z, s.n = z, n
    return false
end

"""
    run_round!(set, chunk) -> Int

Advance every pending sample by `chunk` iterations and swap-remove the ones
that resolved. Returns how many resolved. No sample is advanced twice before
every other has advanced once.
"""
function run_round!(set::SampleSet, chunk)
    resolved = 0
    i = 1
    while i <= length(set.pending)
        idx = set.pending[i]
        if advance!(set.samples[idx], chunk, set.escape_radius2, set.maxiter)
            set.pending[i] = set.pending[end]   # swap-remove
            pop!(set.pending)
            resolved += 1
        else
            i += 1
        end
    end
    return resolved
end

"""
    run_for!(set, budget_ns; chunk = 64) -> Int

Run rounds until `budget_ns` nanoseconds have elapsed or nothing is pending.
This is the whole non-blocking story: work is chopped into slices no longer
than one round, so a caller on the render thread can never stall.
"""
function run_for!(set::SampleSet, budget_ns; chunk = 64)
    resolved = 0
    t0 = time_ns()
    while !isempty(set.pending)
        resolved += run_round!(set, chunk)
        time_ns() - t0 >= budget_ns && break
    end
    return resolved
end

"""
    resolve_all!(set; chunk = 64) -> Int

Run rounds until nothing is pending. For headless use and tests.
"""
function resolve_all!(set::SampleSet; chunk = 64)
    resolved = 0
    while !isempty(set.pending)
        resolved += run_round!(set, chunk)
    end
    return resolved
end

# ---------------------------------------------------------------------------
# Samplers
# ---------------------------------------------------------------------------

"""
    halton(i, base)

The `i`-th element (1-based) of the van der Corput sequence in `base`.
"""
function halton(i, base)
    f = 1.0
    r = 0.0
    while i > 0
        f /= base
        r += f * (i % base)
        i ÷= base
    end
    return r
end

"""
    sprinkle!(set, (xmin, xmax), (ymin, ymax), n; offset = length(set), mindist = 0.0)

Add `n` low-discrepancy (Halton base 2/3) samples over the given rectangle.
Any prefix of the sequence is itself an even cover, so repeated sprinkles
progressively refine instead of clumping. Returns the number actually added.
"""
function sprinkle!(set::SampleSet, xlims, ylims, n;
                   offset = length(set), mindist = 0.0)
    xmin, xmax = xlims
    ymin, ymax = ylims
    added = 0
    for k in 1:n
        i = offset + k
        x = xmin + (xmax - xmin) * halton(i, 2)
        y = ymin + (ymax - ymin) * halton(i, 3)
        push_sample!(set, complex(x, y); mindist) === nothing || (added += 1)
    end
    return added
end

"""
    sprinkle_disc!(set, centre, radius, n; offset = length(set), mindist = 0.0)

Add `n` low-discrepancy samples in the disc of radius `radius` about
`centre`, and return how many were actually added.

`r = radius·√h` rather than `r = radius·h` — without the square root the
samples pile up in the middle, because a disc has more area further out.
"""
function sprinkle_disc!(set::SampleSet, centre, radius, n;
                        offset = length(set), mindist = 0.0)
    c = ComplexF64(centre)
    added = 0
    for k in 1:n
        i = offset + k
        r = radius * sqrt(halton(i, 2))
        θ = 2π * halton(i, 3)
        push_sample!(set, c + r * cis(θ); mindist) === nothing || (added += 1)
    end
    return added
end

"""
    smoothcount(s::Sample, escape_radius2) -> Float64
    smoothcount(set::SampleSet, s::Sample) -> Float64

The renormalised ("smooth") iteration count,

    ν = n + 1 − log₂( log|zₙ| / log R ),     R = √escape_radius2

which is the integer escape count corrected by how far past the escape
radius the orbit actually landed. Returns `NaN` for an `Interior` sample and
the bare `n` for a `Pending` one.

This is §2c.2 of the plan, and it is worth more than de-banding. The integer
count is a staircase: it takes only the values `n`, so *any* interpolating
reconstruction drawn on top of it interpolates a staircase and produces
facets aligned to the escape-time bands. `ν` is continuous across band
boundaries, which is what makes a piecewise-linear reconstruction
([`meshbrotplot!`](@ref)) mean anything at all.

Continuity is worth checking rather than believing. A sample that escapes at
step `n` having barely crossed R gets `ν ≈ n + 1`; one that escapes at `n`
having overshot to `≈ R²` gets `ν ≈ n`. Perturbing that second sample so it
escapes one step earlier, at `n − 1` and barely past R, also gives `ν ≈ n`.
The two branches meet.

`ν ≥ n ≥ 1` for any escaped sample, which keeps the default `log10`
colorscale in its domain.
"""
function smoothcount(s::Sample, escape_radius2)
    s.status === Escaped || return s.status === Interior ? NaN : Float64(s.n)
    az = abs(s.z)
    logR = log(escape_radius2) / 2
    (az > 1 && logR > 0) || return Float64(s.n)
    return s.n + 1 - log(log(az) / logR) / log(2)
end

smoothcount(set::SampleSet, s::Sample) = smoothcount(s, set.escape_radius2)

"""
    value(set, s; smooth = true) -> Float64

The number a sample contributes to the picture: its smooth count, or its raw
integer escape count when `smooth = false`. `NaN` for `Interior` either way.
"""
value(set::SampleSet, s::Sample; smooth = true) =
    smooth ? smoothcount(set, s) :
             (s.status === Interior ? NaN : Float64(s.n))
