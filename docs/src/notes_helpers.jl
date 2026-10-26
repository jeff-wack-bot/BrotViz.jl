# Drawing glue for the figures in notes.md. Trees, layouts, rays, Julia sets
# and wedge diagrams come from BrotViz; this file only adds the small pieces the
# handwritten pages use that BrotViz has no single call for (kneading-word
# labels, the critical partition of the circle, triods and their preimages).

using CairoMakie, BrotViz, Mandelbrot
using Mandelbrot: KneadingSymbol, Sequence

# The A/B colors of BrotViz's wedge diagrams.
const ACOLOR = Makie.to_color(:royalblue)
const BCOLOR = Makie.to_color(:orangered)
const LEGCOLORS = [:red, :green, :blue, :purple]

cis2π(t) = Point2f(cospi(2t), sinpi(2t))

function bare!(ax)
    hidedecorations!(ax)
    hidespines!(ax)
    return ax
end

treeaxis(pos; kwargs...) = bare!(Axis(pos; aspect=DataAspect(), kwargs...))

# Embedded trees live in the dynamical plane; their rays run out to radius 100,
# so pin the view to the Julia set as `embeddedtreeplot` does.
juliaaxis(pos; kwargs...) = treeaxis(pos; limits=(-2, 2, -2, 2), kwargs...)

frac(r::Rational) = "$(numerator(r))/$(denominator(r))"

# An angled internal address as a title: 1_{1/3} → 3_{2/3} → 7.
addresstitle(aia::AngledInternalAddress) = Makie.latexstring(
    join(("$(n)_{$(frac(a))}" for (n, a) in zip(aia.addr, aia.angles)), " \\to ") *
    " \\to $(aia.addr[end])")

# Kneading word of a tree vertex, without the subscript Mandelbrot puts on `*`.
word(node) = join(replace(string(s), r"[₀-₉]" => "") for s in node.items)

# Orbit-index labels, the same ones BrotViz's tree plots use.
orbitlabels(H, nodes) = BrotViz._tree_node_labels(nodes, orbit(H.criticalpoint))

"""
    partialtree(K, n)

The Hubbard tree spanned by the first `n` points of the critical orbit of the
kneading sequence `K`, built by the same triod insertions as `HubbardTree(K)`
but stopped early, so the intermediate trees of the algorithm can be drawn.
"""
function partialtree(K, n)
    K = Sequence{KneadingSymbol{3}}(K)
    starK = Mandelbrot.prepend(K, KneadingSymbol{3}("*"))
    points = orbit(starK).items
    H = Dict(starK => Set([K]), K => Set([starK]))
    for p in points[3:n]
        p in keys(H) || (H = Mandelbrot.addsequence(H, K, (starK, deepcopy(H[starK])), p))
    end
    return HubbardTree(H, starK)
end

"""
    labelledtreeplot!(ax, H; layout, labels=:word)

`treeplot!` with each vertex labelled by its kneading word (`labels=:word`) or
its orbit index (`labels=:orbit`).
"""
function labelledtreeplot!(ax, H; layout=reingold_tilford_layout(H), labels=:word,
                           fontsize=14)
    treeplot!(ax, H, layout; show_labels=false)
    _, nodes = Mandelbrot.adjlist(H.adj)
    pos = [Point2f(p[1], -p[2]) for p in layout.positions]
    text = labels == :word ? word.(nodes) : orbitlabels(H, nodes)
    text!(ax, pos; text, align=(:center, :bottom), offset=(0, 6), fontsize)
    ax.xautolimitmargin = (0.15, 0.15)
    ax.yautolimitmargin = (0.08, 0.15)
    return ax
end

"""
    partitiondisk!(ax, θ; marks=[], labels=("A", "B"))

The circle at infinity cut by the diameter through the two preimages `θ/2`
and `(θ+1)/2` of `θ`. The half containing `θ` gets `labels[1]`, the other
`labels[2]`. `marks` is a list of `angle => label` pairs drawn on the boundary.
"""
function partitiondisk!(ax, θ; marks=Pair[], labels=("A", "B"))
    a, b = θ / 2, (θ + 1) / 2
    poly!(ax, pushfirst!([cis2π(t) for t in range(a, b, 100)], Point2f(0, 0)), color=(ACOLOR, 0.3))
    poly!(ax, pushfirst!([cis2π(t) for t in range(b, a + 1, 100)], Point2f(0, 0)), color=(BCOLOR, 0.3))
    lines!(ax, [cis2π(t) for t in range(0, 1, 200)], color=:black)
    lines!(ax, [cis2π(a), cis2π(b)], color=:black)
    text!(ax, 0.5 * cis2π((a + b) / 2); text=labels[1], align=(:center, :center), fontsize=20)
    text!(ax, 0.5 * cis2π((b + a + 1) / 2); text=labels[2], align=(:center, :center), fontsize=20)
    for (t, label) in marks
        scatter!(ax, [cis2π(t)], color=:black)
        text!(ax, 1.2 * cis2π(t); text=label, align=(:center, :center))
    end
    limits!(ax, -1.45, 1.45, -1.45, 1.45)
    return ax
end

# A leg is a polyline; legs(w, φs) are straight legs of length r leaving w.
legs(w, φs; r=0.35, n=40) = [[w + s * r * cis(φ) for s in range(0, 1, n)] for φ in φs]

# One continuous branch of the pullback of a leg under z ↦ z² + c; the other
# branch is its negative.
function pullback(leg, c)
    z = [sqrt(leg[1] - c)]
    for p in leg[2:end]
        s = sqrt(p - c)
        push!(z, abs(s - z[end]) < abs(s + z[end]) ? s : -s)
    end
    return z
end

function legsplot!(ax, ls; colors=LEGCOLORS, linewidth=2)
    for (leg, col) in zip(ls, colors)
        lines!(ax, real(leg), imag(leg); color=col, linewidth)
    end
    return ax
end

"""
    arms(HC, node; r=0.25)

The arms of the embedded tree `HC` at `node`: the start of each edge leaving
it, oriented outward and cut off at distance `r` from the vertex.
"""
function arms(HC, node; r=0.25)
    z0 = HC.vertices[node]
    out = Vector{ComplexF64}[]
    for nbr in HC.htree.adj[node]
        path = HC.edges[Set([node, nbr])][2]
        abs(path[1] - z0) > abs(path[end] - z0) && (path = reverse(path))
        k = something(findfirst(z -> abs(z - z0) > r, path), length(path))
        push!(out, path[1:k])
    end
    return out
end

"""
    embedded(θ; budget=30_000)
    embedded(aia::AngledInternalAddress; budget=30_000)

The embedded tree that `embeddedtreeplot!(ax, θ)` draws, refined by up to the
same five rounds of `refinetree!` as `HyperbolicComponent(θ)` but stopping once
the edges hold `budget` points. Each round multiplies the edge length by 5 to
30; five rounds of the period-12 tree of 7/13 need over 9 GB.
"""
embedded(θ::Rational; kwargs...) = embedded(OrientedHubbardTree(θ); kwargs...)
embedded(aia::AngledInternalAddress; kwargs...) = embedded(OrientedHubbardTree(aia); kwargs...)

function embedded(OHT::OrientedHubbardTree; budget=30_000)
    HC = HyperbolicComponent(OHT)
    for _ in 1:5
        sum(length(e[2]) for e in values(HC.edges)) > budget && break
        Mandelbrot.refinetree!(HC)
    end
    return HC
end
