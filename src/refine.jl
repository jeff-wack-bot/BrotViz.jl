# Adaptive refinement — letting the data choose where the next samples go.
#
# Every pair of Delaunay-adjacent resolved samples disagrees by some amount.
# Where two neighbours agree, the flat cell between them is already telling
# the truth and another sample there buys nothing; where they disagree, the
# reconstruction is guessing. So rank the Delaunay edges by disagreement and
# drop a new sample at the midpoint of the worst ones.
#
# This is the real generalization of the grid. A grid is a sampling measure
# fixed a priori; this is one the picture asks for. Because ∂M is exactly
# where escape time is discontinuous, samples pile up on it without ever
# being told it exists — and the refinement self-limits, because once two
# neighbours both sit at `maxiter` there is nothing left for a sample
# between them to distinguish.
#
# The triangulation comes from DelaunayTriangulation.jl via the `DelTri`
# alias set up in reconstruction.jl.

"""
    banddisagreement(νa, νb, len)

Default edge priority: the disagreement a *viewer* sees across the edge,
times the length over which they see it.

Two corrections to the bare `|νa − νb|` of the plan's §2a.4, both of which
matter in practice:

  * **Log, not linear.** Escape counts are wildly skewed — over the classic
    view most exterior samples leave in under ten iterations while a few near
    ∂M take hundreds. A linear difference therefore scores the near-boundary
    region orders of magnitude above everything else and spends the entire
    budget there, leaving the escape-time bands in the open exterior as
    coarse as they started. The log difference is also what the default
    `colorscale = log10` actually renders.
  * **Times length.** Disagreement alone is scale-free, so it keeps picking
    the same fractally-disagreeing spot and subdivides it forever. Weighting
    by edge length turns the score into an estimate of the *integrated*
    reconstruction error over the edge, which falls as the edge is halved.
    That is what makes the refinement spread along ∂M instead of boring a
    hole into one point of it.
"""
banddisagreement(νa, νb, len) = abs(log(νa) - log(νb)) * len

"""
    rawdisagreement(νa, νb, len)

§2a.4 as originally written: the bare escape-count difference, ignoring edge
length. Kept because it is the obvious thing to try and it is instructive to
watch it fail — see [`banddisagreement`](@ref).
"""
rawdisagreement(νa, νb, len) = abs(νa - νb)


"""
    planefit_prediction(dx, dy, w) -> Float64

Least-squares fit of `w ≈ α + β·dx + γ·dy` over a vertex's neighbours, with
`dx, dy` measured *from* that vertex, evaluated back at it — so the answer is
just `α`: what linear interpolation of the neighbourhood predicts the vertex
should be. Falls back to the plain mean when the neighbours are collinear (or
too few), which is the best a plane can do there anyway.

Fitting a plane rather than taking a weighted average of neighbours matters:
an average is only exact for a linear field if the vertex sits at the
neighbours' centroid, so on an irregular mesh it would report curvature where
there is only irregularity — and our meshes are deliberately irregular.
"""
function planefit_prediction(dx, dy, w)
    m = length(w)
    m == 0 && return NaN
    s0 = Float64(m)
    sx = sum(dx); sy = sum(dy)
    sxx = sum(x -> x * x, dx); syy = sum(y -> y * y, dy)
    sxy = sum(i -> dx[i] * dy[i], eachindex(dx))
    tw = sum(w); twx = sum(i -> w[i] * dx[i], eachindex(w))
    twy = sum(i -> w[i] * dy[i], eachindex(w))
    # 3x3 symmetric solve by Cramer's rule
    det = s0 * (sxx * syy - sxy * sxy) - sx * (sx * syy - sxy * sy) +
          sy * (sx * sxy - sxx * sy)
    scale = max(s0 * sxx * syy, eps())
    abs(det) <= 1.0e-12 * scale && return tw / s0        # degenerate: use the mean
    return (tw * (sxx * syy - sxy * sxy) - sx * (twx * syy - sxy * twy) +
            sy * (twx * sxy - sxx * twy)) / det
end

"""
    interpolation_residuals(pts, w, isint, neighbours) -> (resid, span)

Per-vertex estimate of how wrong a piecewise-linear reconstruction is *at
that vertex*: the gap between its own value and what a plane through its
Delaunay neighbours predicts. Values are in colour space (`log ν`), because
that is what the mesh actually interpolates — Makie converts per-vertex
scalars to RGBA and interpolates the colours, not the counts.

Three classes, because a discontinuity has no gradient to fit:

  * an interior vertex all of whose neighbours are interior is perfectly
    predicted (residual 0) — the interior is flat and needs almost no samples;
  * a vertex with a neighbour of the other class sits on ∂M, where there is no
    gradient to fit, so it gets the largest *local* jump to a neighbour;
  * everything else gets the plane-fit residual.

Using the local jump rather than the global colour span for that middle case
is not cosmetic. The span is ~500× a typical smooth residual, so boundary
edges outrank every other edge until they are 500× shorter — which starves
the smooth exterior *and* makes refinement fight the pruner, since bisecting
a boundary edge usually drops a point into the flat interior, where the next
pass immediately deletes it again.
"""
function interpolation_residuals(pts, w, isint, neighbours)
    n = length(pts)
    finite = [w[i] for i in 1:n if !isint[i]]
    span = isempty(finite) ? 1.0 : max(maximum(finite) - minimum(finite), eps())
    resid = zeros(n)
    dx = Float64[]; dy = Float64[]; wv = Float64[]
    for v in 1:n
        nb = neighbours[v]
        if any(u -> isint[u] != isint[v], nb)
            resid[v] = maximum(u -> abs(w[v] - w[u]), nb; init = 0.0)
            continue
        end
        isint[v] && continue           # interior among interior: nothing to learn
        empty!(dx); empty!(dy); empty!(wv)
        for u in nb
            push!(dx, pts[u][1] - pts[v][1])
            push!(dy, pts[u][2] - pts[v][2])
            push!(wv, w[u])
        end
        pred = planefit_prediction(dx, dy, wv)
        resid[v] = isfinite(pred) ? abs(w[v] - pred) : span
    end
    return resid, span
end

"""
    refine!(set, n; priority = banddisagreement, mindist = 0.0, within = nothing,
            prune = true, prunetol = 0.0) -> (added = ..., removed = ...)

Move the sample set towards the picture: insert up to `n` new samples at the
midpoints of the most-disagreeing Delaunay edges, and delete the ones that
are not earning their keep.

A sample that the reconstruction could have guessed is redundant: deleting it
does not change the image, so it is not paying for itself. Under
`:disagreement` that means every Delaunay neighbour carries the same value;
under `:interpolation` it means a plane through the neighbours already
predicts it. Deleting every such sample at once is safe, not just one at a
time: a sample next to one the reconstruction gets *wrong* is never
redundant, so the rim of any well-predicted region survives and keeps the
region's outline exactly where it was. This is what
makes the interior cheap — it is 17% of the classic view's area and needs
almost no samples to render correctly.

Addition and deletion cannot fight over the same sample. Redundancy means
*every* incident edge has zero disagreement, while an edge is only worth
bisecting if its disagreement is positive, so the two sets are disjoint by
construction. (That argument is about the raw disagreement, which is what
the redundancy test uses — a custom `priority` that ranks agreeing edges
highly would break the disjointness, though not much else.)

Both halves come out of one triangulation, which is the expensive part.

Only resolved samples carry a value, so only they take part in the
triangulation; pending ones are ignored and will join the next refinement.

`mindist` is forwarded to [`push_sample!`](@ref). **Pass one.** It is not
merely an optimisation: without it a midpoint can land exactly on an
existing sample, DelaunayTriangulation silently drops the duplicate, the
edge that produced it therefore survives into the next triangulation, and
refinement inserts the same dead point on every subsequent call. A tolerance
scaled to the viewport makes that impossible. About two pixels' worth is the
right order — `scatterbrot.jl` passes the same 0.2%-of-width it uses for
clicks — since it also stops refinement from subdividing below what the
screen can show, and it shrinks automatically as you zoom. Rejected midpoints do not consume budget: the walk down the
ranked edges continues until `n` samples have actually been added.

`within = ((xmin, xmax), (ymin, ymax))` restricts refinement to samples in
that rectangle. Worth passing the viewport (padded a little, so that edges
crossing the border are still seen): refining samples the user cannot see is
wasted effort, and the full retriangulation this does on every call is the
cost that grows with the sample count — see the note below.

`priority` is called as `priority(νa, νb, len)`, where `ν` is the iteration
count at which the sample resolved (so an interior sample contributes
`maxiter`) and `len` is the edge length.

`criterion` picks what "not earning its keep" and "worth bisecting" mean:

  * `:disagreement` (default) — the `|Δν|`-based rule, via `priority`.
  * `:interpolation` — how badly a plane through a sample's Delaunay
    neighbours predicts the sample itself
    ([`interpolation_residuals`](@ref)), in colour space. **Use with
    `prune = false`; see the measurements below.**

`:interpolation` is the theoretically-motivated choice for the piecewise-
*linear* reconstruction of [`meshbrotplot!`](@ref): there a region of
constant gradient is reproduced exactly and needs no samples at all, so an
indicator based on how fast the field is *changing* should waste effort in
the smooth exterior, and curvature should be what costs you. Measured
against ground truth on a 250×250 grid, scoring the linear reconstruction
itself, it is not:

| criterion | prune | n | misclassified | linear-err |
|---|---|---|---|---|
| `:disagreement` | yes | 16089 | **0.0025** | **0.0208** |
| `:disagreement` | no  | 16000 | 0.0027 | 0.0175 |
| `:interpolation`| no  | 16000 | 0.0023 | 0.0235 |
| `:interpolation`| yes | 3801  | 0.0049 | 0.0457 |

Two findings, both against the hypothesis. `:interpolation` is marginally
the best at *classifying* the boundary and clearly worse in the smooth
exterior — the opposite of the reason for trying it. And with pruning on it
does not converge: it added 14 000 samples and deleted 10 599 of them to end
at 3 801, because a plane-fit residual measures whether a vertex is
predictable *at the current spacing*, which every vertex becomes as the mesh
refines. So it declares itself finished long before the visible error is
small, while refinement keeps proposing the same edges. `banddisagreement`
happens to be the better indicator for both reconstructions; the honest
reading is that in the exterior the gradient and the curvature of `log ν`
grow together, so the cheaper proxy loses nothing.

`prune = false` turns deletion off. `prunetol` is the residual below which a
sample counts as predictable; `nothing` means 1% of the colour span under
`:interpolation`, and exact equality under `:disagreement`, which is right
for integer escape counts. `smooth` selects raw or renormalised counts for
the `:interpolation` residuals. Samples on the convex hull are never deleted — their cells are the unbounded
ones that reach the clip rectangle, so losing one would visibly redraw the
outer edge, and keeping them also guarantees that what remains can still be
triangulated.

Cost is dominated by retriangulating from scratch: roughly 5 ms at 2 000
resolved samples, 30 ms at 5 000, 190 ms at 20 000 and over a second at
60 000 — superlinear, because samples clustered on a fractal curve are close
to the worst case for Bowyer–Watson (letting it randomise the insertion
order barely helps; it was measured). At interactive sample counts this is
imperceptible; past ~20 000 it is a visible hitch on the keypress. The fix is
to maintain the triangulation incrementally with `DelTri.add_point!` instead
of rebuilding it, which §9 Phase 4 already owns.
"""
function refine!(set::SampleSet, n; criterion = :disagreement,
                 priority = banddisagreement, mindist = 0.0, within = nothing,
                 prune = true, prunetol = nothing, smooth = true)
    resolved = Int[]
    for i in eachindex(set.samples)
        s = set.samples[i]
        s.status === Pending && continue
        if within !== nothing
            (xmin, xmax), (ymin, ymax) = within
            xmin <= real(s.c) <= xmax && ymin <= imag(s.c) <= ymax || continue
        end
        push!(resolved, i)
    end
    length(resolved) < 3 && return (added = 0, removed = 0)

    pts = [(real(set.samples[i].c), imag(set.samples[i].c)) for i in resolved]
    triangulable(pts) || return (added = 0, removed = 0)

    tri = DelTri.triangulate(pts; randomise = false)

    neighbours = [Int[] for _ in eachindex(pts)]
    edges = Tuple{Int, Int}[]
    for (a, b) in DelTri.each_solid_edge(tri)
        push!(neighbours[a], b)
        push!(neighbours[b], a)
        push!(edges, (a, b))
    end

    scores = Float64[]
    midpoints = ComplexF64[]
    keepvertex = falses(length(pts))     # true = not redundant

    if criterion === :interpolation
        isint = [set.samples[i].status === Interior for i in resolved]
        # An interior sample gets log(maxiter): not a real count, but the right
        # scale, so that the jump across ∂M is a comparable number rather than
        # an infinity.
        capped = log(max(set.maxiter, 2))
        w = [isint[k] ? capped :
             log(max(value(set, set.samples[resolved[k]]; smooth), 1.0))
             for k in eachindex(resolved)]
        resid, span = interpolation_residuals(pts, w, isint, neighbours)
        tol = prunetol === nothing ? 0.01 * span : prunetol
        for v in eachindex(pts)
            keepvertex[v] = resid[v] > tol
        end
        for (a, b) in edges
            pa, pb = pts[a], pts[b]
            len = hypot(pb[1] - pa[1], pb[2] - pa[2])
            len > 0 || continue
            # error estimate on the edge, times the length it spans
            push!(scores, max(resid[a], resid[b]) * len)
            push!(midpoints, complex((pa[1] + pb[1]) / 2, (pa[2] + pb[2]) / 2))
        end
    else
        tol = prunetol === nothing ? 0.0 : prunetol
        ν = [Float64(set.samples[i].n) for i in resolved]
        spread = zeros(length(pts))     # largest disagreement each sample has with a neighbour
        for (a, b) in edges
            d = abs(ν[a] - ν[b])
            spread[a] = max(spread[a], d)
            spread[b] = max(spread[b], d)
            pa, pb = pts[a], pts[b]
            len = hypot(pb[1] - pa[1], pb[2] - pa[2])
            # A zero-length edge's midpoint *is* an existing sample; inserting it
            # would be dropped as a duplicate and the edge would come back
            # forever. `mindist` covers this too, but not everyone passes one.
            len > 0 || continue
            push!(scores, priority(ν[a], ν[b], len))
            push!(midpoints, complex((pa[1] + pb[1]) / 2, (pa[2] + pb[2]) / 2))
        end
        for v in eachindex(pts)
            keepvertex[v] = spread[v] > tol
        end
    end

    # Delete before inserting, so the doomed indices are still valid.
    removed = 0
    if prune
        hull = Set(DelTri.get_vertices(DelTri.get_convex_hull(tri)))
        doomed = [resolved[v] for v in eachindex(pts)
                  if !keepvertex[v] && !(v in hull)]
        if !isempty(doomed)
            deleteat!(set.samples, doomed)      # ascending, as `deleteat!` needs
            # `pending` holds indices, and every one past a deletion has moved
            empty!(set.pending)
            for i in eachindex(set.samples)
                set.samples[i].status === Pending && push!(set.pending, i)
            end
            removed = length(doomed)
        end
    end

    added = 0
    for k in sortperm(scores; rev = true)
        added >= n && break             # checked first, so `n = 0` adds nothing
        scores[k] > 0 || break          # the rest agree exactly; nothing to learn
        push_sample!(set, midpoints[k]; mindist) === nothing || (added += 1)
    end
    return (added = added, removed = removed)
end
