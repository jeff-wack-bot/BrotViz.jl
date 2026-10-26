# Mandelbrot and BrotViz development notes

August 2, 2026

```@setup notes
include("notes_helpers.jl")
```

## The four embeddings of the internal address [1,3,7]

An example: ``[1,3,7]`` has 4 embeddings which correspond to placing vertex ``1`` on any of these 4 vertices where are ``\beta`` and ``-\beta``?

```math
\begin{aligned}
0.001010* &\;\rightsquigarrow\; 0.\overline{0010011} && 19/127\\
&\phantom{\;\rightsquigarrow\;} 0.\overline{0010100} && 20/127\\
0.010001* &\;\rightsquigarrow\; 0.\overline{0100011} && 35/127\\
&\phantom{\;\rightsquigarrow\;} 0.\overline{0100100} && 36/127
\end{aligned}
```

```math
0.010100* \;\; ? \qquad \text{(crossed out on the page)}
```

```math
0.100010* \;\; ? \qquad \text{(crossed out on the page)}
```

Does changing the embedding change the location of ``\beta`` and ``-\beta``?

```math
\begin{aligned}
1_{1/3} \to 3_{1/3} \to 7 && 19/127 \ \&\ 20/127 && 0010011,\; 0010100\\
1_{1/3} \to 3_{2/3} \to 7 && 35/127 \ \&\ 36/127 && 0100011,\; 0100100\\
1_{2/3} \to 3_{1/3} \to 7 && 91/127 \ \&\ 92/127 && 1011011,\; 1011100\\
1_{2/3} \to 3_{2/3} \to 7 && 107/127 \ \&\ 108/127 && 1101011,\; 1101100
\end{aligned}
```

Written loose on the sheet, to the right of the first pair of rows: ``010001*``. Written under the third tree: ``001100*``, ``110011*``.

*The only two good embeddings?* *(green note beside the first pair of trees)*

*Somehow this/these embedding is bad because the branch point dynamics are wrong! characteristic point needs AAB* *(red note around the second tree)*

*This is wrong too? Why?* *(red note under the third tree)*

```@example notes
# The four angled internal addresses over [1,3,7], two per figure. Each is
# drawn from its explicit address: `AngledInternalAddress(θ)` returns 1/3 at
# every level, so `embeddedtreeplot!(ax, 35//127)` would draw 19/127's tree.
fig = Figure(size=(800, 400))
for (j, angles) in enumerate(([1//3, 1//3], [1//3, 2//3]))
    aia = AngledInternalAddress([1, 3, 7], angles)
    embeddedtreeplot!(juliaaxis(fig[1, j], title=addresstitle(aia)), aia)
end
fig
```

*Two Hubbard-tree embeddings for the internal address (1,3,7), annotated in green as the only two good ones*

```@example notes
fig = Figure(size=(450, 400))
aia = AngledInternalAddress([1, 3, 7], [2//3, 1//3])
embeddedtreeplot!(juliaaxis(fig[1, 1], title=addresstitle(aia)), aia)
fig
```

*An embedding of the (1,3,7) tree annotated in red as bad because the branch point dynamics are wrong*

```@example notes
fig = Figure(size=(450, 400))
aia = AngledInternalAddress([1, 3, 7], [2//3, 2//3])
embeddedtreeplot!(juliaaxis(fig[1, 1], title=addresstitle(aia)), aia)
fig
```

*A further embedding of the tree, annotated in red as wrong too*

## The 1/11 Hubbard tree and why the current algorithm works on it

Example ``1/11``

``15`` is a preimage of ``13``, it has fewer arms and nothing will map to the ``2`` arm.

```math
\begin{array}{ccccccccccl}
1&2&3&4&5&6&7&8&9&10&\\
0&0&0&1&0&1&1&1&1&0& {}=94/1023\\
0&0&0&1&0&1&1&1&0&1& {}=1/11
\end{array}
```

Current alg works on ``1/11`` because preimage

*A red curve is drawn through the tree separating a region marked ``0`` from a region marked ``1``.*

```@example notes
# `embedded` (in notes_helpers.jl) is HyperbolicComponent(1//11) with the
# refinement capped; see its docstring.
fig = Figure(size=(600, 600))
embeddedtreeplot!(juliaaxis(fig[1, 1], title="1/11"), embedded(1//11))
fig
```

*Hubbard tree for the external angle 1/11, vertices numbered 0 to 18 with beta and -beta endpoints*

## Building the 1/11 tree: characteristic points, edge list, cyclic order

Example: ``1/11 = 0001011101``

```math
\begin{array}{l}
0001110101 = \\
0011100101\\
\qquad\qquad\quad 10
\end{array}
```

The tree is marked ``A`` above the spine and ``B`` below it.

We have an edge list

```math
\begin{array}{c|l}
\text{vertex} & \text{neighbors}\\\hline
1 & 11\\
2 & 11\\
3 & 12\\
4 & 8\\
5 & 12\\
6 & 14\\
7 & 14\\
8 & 4,\,13
\end{array}
\qquad\qquad
\begin{array}{c|l}
9 & 11,\,12\\
10 & 11,\,13\\
11 & 1,\,2,\,9,\,10\\
12 & 3,\,5,\,9\\
13 & 8,\,10,\,14\\
14 & 6,\,7,\,13
\end{array}
```

1. Identify the characteristic points
2. use the numerators to put a cyclic order on the arms of the charachteristic points
3. use the dynamics to put cyclic order on the other branch points.

add a cyclic order list, same shape as edge list

```@example notes
# Vertices off the critical orbit are labelled by their kneading sequences.
fig = Figure(size=(800, 400))
treeplot!(treeaxis(fig[1, 1], title="1/11"), HubbardTree(KneadingSequence(1//11)); layout=:rt)
fig
```

*The 1/11 tree redrawn flat on a horizontal spine, with kneading words ABA, BA, BBA at the branch vertices*

## Two embeddings of the 7/13 tree

Example: ``7/13``

```math
001001110100
```

First embedding: a horizontal spine from ``-\beta`` on the left to ``\beta`` on the right, carrying the marked points ``20,\ 19,\ 7,\ 11,\ 15,\ 12=0,\ 18,\ 22,\ 21,\ 23``, with arms to ``5,\ 9,\ 2`` above and ``3,\ 10,\ 8,\ 7,\ 6`` below, and a subtree ``4 - 16 - 1`` rising above the point ``15``.

```math
010001110110
```

Second embedding: the same spine with marked points ``20,\ 19,\ 17,\ 11,\ 15,\ 12=0,\ 18,\ 22,\ 21,\ 23``, but the arm at the third point now hangs below to ``2``, and the subtree above ``15`` is ``1 - 16 - 4``, i.e. the two leaves are exchanged.

Beside the leaf ``1`` of the second subtree: “after 1”.

```@example notes
fig = Figure(size=(600, 600))
embeddedtreeplot!(juliaaxis(fig[1, 1], title="7/13"), embedded(7//13))
fig
```

*First embedding of the 7/13 tree drawn on a horizontal spine from -beta to beta*

```@example notes
# 7/13 has angled internal address 1 → 2 → 3 → 5 → 6 → 12 with every angle
# 1/2, so its tree has only the one embedding above. Shown here is the
# abstract tree, before any cyclic order is chosen.
fig = Figure(size=(800, 400))
treeplot!(treeaxis(fig[1, 1], title="7/13, unembedded"), HubbardTree(KneadingSequence(7//13)); layout=:rt)
fig
```

*Second embedding of the 7/13 tree, with the subtree leaves and the third arm exchanged*

## Orbit pictures for 5/11 and 4/9, and a Hubbard tree with dynamics

Another example: ``5/11``

```math
01110\,\overbrace{0}^{6}\,1110
```

the area “flipped” seems to be between ``-\beta`` and ``\alpha``

```math
001100
```

``4/9``

```math
011100
```

```@example notes
fig = Figure(size=(600, 450))
embeddedtreeplot!(juliaaxis(fig[1, 1], title="5/11"), 5//11)
fig
```

*The orbit of 5/11 laid out along an interval from ``-\beta`` to ``\beta``*

```@example notes
fig = Figure(size=(600, 450))
embeddedtreeplot!(juliaaxis(fig[1, 1], title="4/9"), 4//9)
fig
```

*The orbit of 4/9 laid out along an interval from ``+\beta`` to ``-\beta``*

```@example notes
# Vertices are coloured by A/B and by side of the (0, 1) boundary.
fig = Figure(size=(600, 450))
treeplot!(treeaxis(fig[1, 1], title="4/9"), 4//9)
fig
```

*Tree of straight segments with numbered vertices and coloured arrows showing the dynamics*

## Matching the two star points to the half-planes; what a Hubbard tree data structure needs

we u

```math
*_1,\; *_2 \;=\; \pm\sqrt{\texttt{Legs[1][End]}/\lambda}
```

we know ``+`` is in RHP and ``-`` is in LHP

what we want to do is match up

```math
\begin{array}{cc}
*_1 & *_+ \\[2pt]
*_2 & *_-
\end{array}
```

when is ``*_1`` in the RHP?

---

we want a more verbose/detailed description of a Hubbard tree. In addition to an adjacency matrix, we want a description of the dynamics. The dynamics are implicit in the kneading sequences of each vertex and are given by the shift map. How do we represent the dynamics in easily readable form?

Question: are the dynamics of branch points ever not a sequence? Yes. ``1/11``

Is ther an example with a smaller tree/shorter kneading sequence? can't find one for now

what data structure can do both?

Its a digraph, but every node has only one arrow coming out of it.

## Example 14/93: characteristic points and the preimage walk

Example: ``14/93``

```math
0010100001
```

Preimage of 14: 16

what are we attatching to 16? 3

so that when ``\beta`` stays and ``0`` goes to ``1``, ``3`` is after ``0`` and before ``\beta``

Preimage of 16: 15

what are we attatching to 15? 2

so when ``\alpha`` stays and ``-\beta`` goes to ``\beta``

done!

---

A characteristic point separates the critical value from ``0`` and from all the forward images of the critical point.

1. Identify the periodic cycles of branch points, excluding ``\beta \to \beta``)
2. for each cycle, identify the characteristic point

   ```math
   \operatorname{iteratetriod}(z_c,\, z,\, 1) \;=\; \text{``flat''}
   ```

   for all ``z`` in the forward orbit.

```@example notes
fig = Figure(size=(600, 600))
embeddedtreeplot!(juliaaxis(fig[1, 1], title="14/93"), embedded(14//93))
fig
```

*Hubbard tree for 14/93 drawn along a red spine from ``-\beta`` to ``\beta``*

## Adding a shift vector to the Hubbard tree; example 91/127 and death to the conjecture

To the Hubbard tree ``(A, [\kappa])`` we add a vector specifying the image of each vertex under the shift map. Done

---

Plotting a tree with cyclic order info:

1. how many total vertices in this generation?

---

Example:

```math
\left.\begin{array}{l} 91/127 \\ 92/127 \end{array}\right\} \;=\; 1011100
```

*start with critical value. vertex 1.*

Death to the conjecture.

```math
\begin{aligned}
3/19 &\to \text{2 non-2 denominators} \\
     &\phantom{{}\to{}} \text{3 branch points between 0 \& 1} \\
5/19 &\to \text{3 non-2 denominators} \qquad \text{(crossed out on the page)} \\
7/19 &\to \text{1 non-2 denominator} \\
     &\phantom{{}\to{}} \text{4 branch points between 0 \& 1} \qquad \text{and no } \alpha!
\end{aligned}
```

```@example notes
fig = Figure(size=(600, 600))
aia = AngledInternalAddress([1, 3, 7], [2//3, 1//3])   # 91/127 and 92/127
embeddedtreeplot!(juliaaxis(fig[1, 1], title="91/127"), aia)
fig
```

*Hubbard tree for 91/127 drawn on a red spine with branches above and below*

```@example notes
# The vertices of the 91/127 tree projected onto the real axis.
HC = embedded(aia)
nodes = collect(keys(HC.vertices))
x = [real(HC.vertices[n]) for n in nodes]
fig = Figure(size=(800, 160))
ax = bare!(Axis(fig[1, 1]))
lines!(ax, [extrema(x)...], [0, 0], color=:red)
scatter!(ax, x, zero(x), color=:black)
# Alternate labels above and below so neighbours do not collide.
order = sortperm(x)
for (k, i) in enumerate(order)
    text!(ax, x[i], 0; text=orbitlabels(HC.htree, nodes)[i], fontsize=11,
        offset=(0, isodd(k) ? 6 : -6), align=(:center, isodd(k) ? :bottom : :top))
end
fig
```

*The same vertices collapsed onto a single red segment*

## Preserving cyclic order of spider legs; regions A and B at infinity

**Fact:** the spider map preserves the cyclic order of the legs. We want to use this fact to cook up a method for selecting image legs which preserves their cyclic order by default.

**Motivation:** when watching gifs of runs that fail to converge legs are seen to be “hopping around.”

**Idea:** store a ‘body’ of the spider in memory. This consists of a disk at infinity.

**Question:** do we use ``\boxed{\texttt{Legs[i][end]}}`` or ``\texttt{newlegs[i][end]}`` to define regions ``A`` and ``B`` at infinity? More generally, the preimages of which of these defines the boundary?

``\underline{\text{No mixing generations!}}``

if `cregion` ``= A`` then

```math
*_2 = \frac{a}{2}, \qquad *_1 = \frac{a+1}{2}
```

else (`cregion` ``= B``)

```math
*_1 = \frac{a}{2}, \qquad *_2 = \frac{a+1}{2}
```

Only one of the halves will be in the bottom half plane, and it will have a greater angle.

```@example notes
# cregion = A: the critical value θ₁ lies in A, and *₂ = a/2, *₁ = (a+1)/2.
a = 1//7
fig = Figure(size=(450, 450))
partitiondisk!(treeaxis(fig[1, 1], title="cregion = A"), a;
    marks=[a => "θ₁", a/2 => "*₂", (a+1)/2 => "*₁", 0 => "0"])
fig
```

*Case cregion = A: two labelled discs with a chord separating regions A and B, and the resulting formulae for the starred angles*

```@example notes
# cregion = B: A and B trade places, and so do *₁ and *₂.
fig = Figure(size=(450, 450))
ax = treeaxis(fig[1, 1], title="cregion = B")
partitiondisk!(ax, a; labels=("B", "A"),
    marks=[a => "θ₁", a/2 => "*₁", (a+1)/2 => "*₂", 0 => "0"])
fig
```

*Case cregion = B: the same pair of discs with A and B interchanged, and the swapped formulae for the starred angles*

## Are all branch points between 0 and 1 characteristic? A simple example

**Question:** are all branch points between ``0`` and ``1`` characteristic?

---

A simple example. It would be very nice to be thinking in sequences.

Focus on ordering all the sequences based on the ordering of the characteristic points.

It should be done with preimages.

Preimages of characteristic points account for all the branch points, forward images of characteristic points do not.

```@example notes
fig = Figure(size=(700, 400))
treeplot!(treeaxis(fig[1, 1], title="1/11"), HubbardTree(KneadingSequence(1//11)); layout=:rt)
fig
```

*A small tree with two branch points, illustrating characteristic points and their preimages*

## Attaching beta and minus beta via the triod map; representing the embedded tree

**Idea for implementing 1):** we know the kneading sequences of ``\beta`` and ``-\beta``. Can we use the triod map to attatch them to the tree? This works! (I think)

---

**Question:** Are all branch points between ``0`` & ``1`` characteristic? Seems to be true… no counterexamples found yet. Nope. ``19/43`` (crossed out on the page) *not a counterexample*

---

How do we represent the embedded tree? To each branch point we attatch an ordering to the arms.

**Idea:** the arm towards zero is always in position ``0 =`` degree of branch; the numerator specifies the position of the arm towards ``1``. The fact that the first return map sends the arm towards ``0`` to the arm towards ``1`` and permutes the arm transitively allows us to order the arms on the characteristic orbit.

```@example notes
fig = Figure(size=(600, 600))
embeddedtreeplot!(juliaaxis(fig[1, 1], title="1/11"), embedded(1//11); show_rays=false)
fig
```

*An embedded tree with orderings attached at each branch point, with a red curve traced through it*

## Embedding a Hubbard tree: cyclic orders at branch points and the current data representation

To embed a Hubbard tree we have to give a cyclic order to each branch point.

**Question:** Are there Hubbard trees with disjoint orbits of branch points? Yes. ``19/127``

Given a Hubbard tree, current data representation

```math
(\text{Adjacency matrix}, \ \text{Kneading sequence list})
```

How do we find the orbit structure of the vertices?

Two things to focus on:

1. Take a Hubbard tree and append ``\beta`` and ``-\beta``.
2. Take a list of numerators and generate a cyclic order on the endpoints.

Leads to

``13.6/13.9`` -- for each point on the critical orbit, what side of the line ``(-\beta, \beta)`` is it on?

## Hubbard trees and the majority-vote chop procedure

What's the deal with Hubbard trees?

A tree and a map which is surjective. The critical orbit has itinerary ``\overline{*\nu}``.

example: ``\nu = \overline{10110*}``

```math
\begin{array}{ccc}
\overline{*10110} & \overline{10110*} & \overline{0110*1}\\
\nu_1 & \nu_2 & \nu_3\\[6pt]
\overline{110*10} & \overline{10*101} & \overline{0*1011}\\
\nu_4 & \nu_5 & \nu_6
\end{array}
```

*What is ``*BACAAB`` ??*

Majority vote ``\to``

```math
\begin{array}{rcccccc}
\nu_1 = A & * & 1 & 0 & 1 & 1 & 0\\
\nu_2 = B & 1 & 0 & 1 & 1 & 0 & *\\
\nu_3 = C & 0 & 1 & 1 & 0 & * & 1
\end{array}
\qquad \text{no chop}\qquad \text{FLAT $\nu_1$ in the middle}
```

```math
\begin{array}{crcccccc}
 & \nu_2 = A & 1 & 0 & 1 & 1 & 0 & *\\
1 & \nu_3 = B & 0 & 1 & 1 & 0 & * & 1\\
 & \nu_4 = C & 1 & 1 & 0 & * & 1 & 0
\end{array}
\qquad \text{chop $B$}
```

```math
\begin{array}{ccccccc}
 & 0 & 1 & 1 & 0 & * & 1\\
1 & 1 & 0 & 1 & 1 & 0 & *\\
 & 1 & 0 & * & 1 & 0 & 1
\end{array}
\qquad \text{chop $A$}
```

```math
\begin{array}{ccccccc}
 & 1 & 0 & 1 & 1 & 0 & *\\
0 & 0 & 1 & 1 & 0 & * & 1\\
 & 0 & * & 1 & 0 & 1 & 1
\end{array}
```

(continuing in the right-hand column)

```math
\begin{array}{ccccccc}
 & 1 & 0 & 1 & 1 & 0 & *\\
1 & 1 & 1 & 0 & * & 1 & 0\\
 & * & 1 & 0 & 1 & 1 & 0
\end{array}
\qquad \text{chop $C$}
```

```math
\begin{array}{ccccccc}
 & 0 & 1 & 1 & 0 & * & 1\\
1 & 1 & 0 & * & 1 & 0 & 1\\
 & 1 & 0 & 1 & 1 & 0 & *
\end{array}
\qquad \text{chop $A$}
```

```math
\begin{array}{ccccccc}
 & 1 & 0 & 1 & 1 & 0 & *\\
0 & 0 & * & 1 & 0 & 1 & 1\\
 & 0 & 1 & 1 & 0 & * & 1
\end{array}
\qquad \text{chop $A$}
```

```math
\begin{array}{ccccccc}
 & 1 & 0 & 1 & 1 & 0 & *\\
1 & * & 1 & 0 & 1 & 1 & 0\\
 & 1 & 1 & 0 & * & 1 & 0
\end{array}
\qquad \text{chop $B$}
```

We conclude ``(1,2,3)`` is flat, 1 in middle; ``(2,3,4)`` is branched.

## Adjacency from flat triples; the map sends branch points to branch points

```math
\begin{array}{cccccc}
0 & 1 & 1 & 0 & * & 1\\
1 & 0 & 1 & 1 & 0 & *\\
1 & 0 & * & 1 & 0 & 1
\end{array}
```

Every triple of sequences is tested (can we get away with less?)

``\ldots`` (crossed out on the page)

If there exists a sequence ``b`` (“between”) such that ``(x,y,b)`` is flat with ``b`` in the middle, i.e. it reaches the stop case with ``b`` “on ``*``”, then ``x`` and ``y`` are not adjacent.

``(1,2,3)`` is flat with 1 in the middle

```math
\Rightarrow \quad \text{2 is not adjacent to 3}
```

``(2,3,4)`` is branched

- 2 is not adjacent to 3
- 3 is not adjacent to 4
- 4 is not adjacent to 2

The map sends branch points to branch points.

```@example notes
# The critical cycle (black) and the cycle of branch points (red) of 1/11.
H = HubbardTree(KneadingSequence(1//11))
layout = reingold_tilford_layout(H)
fig = Figure(size=(700, 400))
ax = treeaxis(fig[1, 1], title="1/11")
labelledtreeplot!(ax, H; layout, labels=:orbit)
E, nodes = Mandelbrot.adjlist(H.adj)
branch = [i for i in eachindex(nodes) if length(E[i]) >= 3]
scatter!(ax, [Point2f(layout.positions[i][1], -layout.positions[i][2]) for i in branch],
    color=:red, markersize=14)
fig
```

*Stacked bracket schematic: critical cycle above branch-points necklace, with a second empty bracket alongside*

## Triod algorithm run on ``*ABABAAA``: regions A and B empty

```math
\begin{aligned}
\text{start:} \qquad & *ABABAAA
\end{aligned}
```

get: ``ABABAAA*``

Region ``A`` is empty.

result:

```math
\begin{array}{c}
ABABAAA*\\
|\\
*ABABAAA
\end{array}
```

get: ``BABAAA*A``

Region ``B`` is empty.

result:

```math
\begin{array}{c}
ABABAAA*\\
|\\
*ABABAAA\\
|\\
BABAAA*A
\end{array}
```

## Running triod``(ABAAA*AB, ABABAAA*, *ABABAAA)``: branched ``\overline{AB}``

get: ``ABAAA*AB``

run triod``(ABAAA*AB,\; ABABAAA*,\; *ABABAAA)``

Majority:

```math
\begin{array}{lll}
\text{BAAA*ABA} & \text{BABAAA*A} & \text{ABABAAA*}\\
\text{AAA*ABAB} & \text{ABAAA*AB} & \text{ABABAAA*}\\
\text{AA*ABABA} & \text{BAAA*ABA} & \text{BABAAA*A}\\
\text{ABABAAA*} & \text{AAA*ABAB} & \text{ABAAA*AB}
\end{array}
\qquad
\begin{array}{c}
A\\ B\\ A\\ B\\ A\\ B\\ A\\ B
\end{array}
```

branched: ``\overline{AB}``

result:

```@example notes
# The triod algorithm for K = ABABAAA* (the kneading sequence of 101/255),
# stopped after inserting ABAAA*AB: the first branch point, AB, appears.
K = KneadingSequence(101//255)
fig = Figure(size=(600, 350))
labelledtreeplot!(treeaxis(fig[1, 1]), partialtree(K, 4))
fig
```

*Resulting triod: two kneading words meeting a third at the branch point ``AB``*

## Running triod on *ABABAAA: majority column and result trees

get: ``BAAA*ABA``

```math
\text{run triod}(\,*ABABAAA \quad BABAAA*A \quad BAAA*ABA\,)
```

```math
\begin{array}{r|ccc}
\text{majority} & & & \\
B & \text{(votes on the \texttt{run triod} line above)} & & \\
A & ABABAAA* & ABAAA*AB & AAA*ABAB\\
B & BABAAA*A & BAAA*ABA & AA*ABABA\\
A & ABAAA*AB & AAA*ABAB & ABABAAA*
\end{array}
```

result:

get: ``AAA*ABAB`` result:

run triod

```math
\begin{array}{cccccccc}
* & A & B & A & B & A & A & A\\
A & B & A & B & A & B & & \\
A & A & A & * & A & B & A & B
\end{array}
```

```@example notes
# ... after also inserting BAAA*ABA: the second branch point, BA.
fig = Figure(size=(700, 400))
labelledtreeplot!(treeaxis(fig[1, 1]), partialtree(K, 5))
fig
```

*Result tree: two sequences at the top meeting at AB, descending through *ABABAAA to BA and branching to two sequences*

```@example notes
# ... and after AAA*ABAB.
fig = Figure(size=(700, 450))
labelledtreeplot!(treeaxis(fig[1, 1]), partialtree(K, 6))
fig
```

*Compact redrawing of the result tree at the bottom right of the page*

```@example notes
# The three sequences of the second triod, column by column.
rows = ["*ABABAAA", "ABABAB", "AAA*ABAB"]
fig = Figure(size=(500, 200))
ax = bare!(Axis(fig[1, 1], aspect=DataAspect()))
for (i, r) in enumerate(rows), (j, ch) in enumerate(r)
    text!(ax, j, -i; text=string(ch), align=(:center, :center), fontsize=i == 2 ? 14 : 24)
end
for j in 1:4
    poly!(ax, Rect2f(j - 0.45, -3.5, 0.9, 3), color=:transparent, strokecolor=:gray50, strokewidth=1)
end
limits!(ax, 0, 9, -4, 0)
fig
```

*Boxed column-by-column comparison used to run the second triod*

## Shower perspective: mapping triods with sqrt

shower perspective: think about mapping triods with sqrt

(The remainder of the page is three drawn triods; see figures.)

```@example notes
# A triod in the dynamical plane of the rabbit, z ↦ z² + c ...
c = parameter(1//7, 250)
w = 0.4 + 0.5im
triod = legs(w, (3π/4, π/4, -π/2))
fig = Figure(size=(400, 400))
ax = treeaxis(fig[1, 1])
legsplot!(ax, triod)
fig
```

*Upright triod drawn with three coloured legs meeting at a central vertex*

```@example notes
# ... and its two preimages ±√(z - c), on either side of the critical point 0.
fig = Figure(size=(400, 400))
ax = treeaxis(fig[1, 1])
legsplot!(ax, [pullback(leg, c) for leg in triod])
legsplot!(ax, [-pullback(leg, c) for leg in triod])
scatter!(ax, [0], [0], color=:black, markersize=12)
fig
```

*Two rotated triods separated by a heavy dot in the lower right of the page*

## Tree traversal and the before/after flip for 2452/16383

Left ``\to 0``, Right ``\to 1`` top ``\to 0``, bottom ``\to 1`` Left ``\to 1``, Right ``\to 0``

The transition from ``7`` to ``8`` flips us from “before” to after. Why?

```math
0010011 \;\big|\; 0010\,\underbrace{100}_{\phantom{x}} \;=\; \frac{2452}{16383}
```

The flip happens when parent goes to child and child goes to parent.

``7`` is where the flip occurs, but it doesn't take effect until we reach ``12``, a non-leaf.

```@example notes
fig = Figure(size=(800, 450))
treeplot!(treeaxis(fig[1, 1], title="2452/16383"), HubbardTree(KneadingSequence(2452//16383)))
fig
```

*Numbered binary tree with green and red traversal paths overlaid, together with the left/right/top/bottom bit legend*

## A hard one: the angle 1/55

Lets do a hard one: ``1/55``

*focus on* (see the star figure)

*no maybes, coprime and smaller*

Answer: ``2``!

*no periodic branch points, only preimages of ``\alpha``*

```@example notes
fig = Figure(size=(800, 600))
treeplot!(treeaxis(fig[1, 1], title="1/55"), HubbardTree(KneadingSequence(1//55)); layout=:rt)
fig
```

*Numbered tree for the angle 1/55, with the alpha fixed point marked and red image-labels on the interior vertices*

```@example notes
# 1/63 = 1/(2⁶ - 1): the star with six spokes.
fig = Figure(size=(450, 450))
embeddedtreeplot!(juliaaxis(fig[1, 1], title="1/63"), 1//63)
fig
```

*Star of six spokes at a single vertex, labelled 0 through 5, under the red heading focus on*

## Vertex-labelled trees for the angles 3/7 and 4/9

Two vertex-labelled trees drawn on a ruled index card, with no accompanying prose. Each carries a green identification of the angle as a periodic binary expansion.

The first tree is annotated

```math
0.\overline{011} = \frac{3}{7},
```

and has eight marked points, labelled ``0,1,2,3,4,5,6,7``, with ``0`` at the apex-adjacent centre and a dashed vertical through it.

The second tree is annotated

```math
0.\overline{011100} = \frac{4}{9},
```

and has ten marked points, labelled ``1,\ldots,10``.

```@example notes
fig = Figure(size=(450, 450))
embeddedtreeplot!(juliaaxis(fig[1, 1], title="0.(011) = 3/7"), 3//7)
fig
```

*Tree for the angle 3/7, vertices labelled 0-7, with the green annotation ``0.\overline{011}=3/7``.*

```@example notes
fig = Figure(size=(450, 450))
embeddedtreeplot!(juliaaxis(fig[1, 1], title="0.(011100) = 4/9"), 4//9)
fig
```

*Tree for the angle 4/9, vertices labelled 1-10, with the green annotation ``0.\overline{011100}=4/9``.*

## A twelve-point tree drawn as a real spine and as a planar tree

The same twelve marked points are drawn twice: first strung out along a single horizontal segment, then re-drawn as a planar tree with a branch point.

Along the horizontal segment the points appear in the order

```math
10,\; 1,\; 5,\; 7,\; 3,\; 12,\; 4,\; 8,\; 6,\; 2,\; 11,\; 9,
```

with the endpoints marked ``-\beta`` (at ``10``) and ``\beta`` (at ``9``), and the interior points ``12`` and ``11`` marked ``\alpha`` and ``-\alpha`` respectively.

Between the two drawings stands the binary block

```math
01101001 .
```

In the second drawing the branch point is ``8``; the left branch runs ``4,\,12\,(\alpha),\,3,\,7,\,5,\,1,\,10`` and the right branch runs ``6,\,2,\,11\,(-\alpha),\,9\,(\beta)``.

```@example notes
# 0.(01101001) = 105/255 = 7/17
fig = Figure(size=(600, 450))
embeddedtreeplot!(juliaaxis(fig[1, 1], title="7/17"), 7//17)
fig
```

*The twelve marked points arranged along a single horizontal spine, from ``-\beta`` on the left to ``\beta`` on the right.*

```@example notes
fig = Figure(size=(600, 450))
treeplot!(treeaxis(fig[1, 1], title="7/17"), 7//17)
fig
```

*The same twelve points redrawn as a planar tree branching at the point 8.*

## Hubbard tree for the internal address 1 -> 2 -> 6 -> 11

Hubbard tree for

```math
1 \to 2 \to 6 \to 11
```

with the angles ``\tfrac{1}{2}``, ``\tfrac{1}{3}``, ``\tfrac{1}{3}`` written beneath the three arrows, the second and third each enclosed in a box.

The kneading sequence and the corresponding angle are recorded as

```math
\overline{ABABAAABABB}, \qquad
0.\overline{01011001011} = \frac{715}{2047}.
```

The tree itself is drawn across the rest of the sheet. Its eleven orbit points are labelled ``1,\ldots,11`` together with the endpoints ``\beta`` and ``-\beta`` and the critical point ``0``; the branches carry the periodic labels

```math
\overline{AB},\quad \overline{ABABAA},\quad \overline{ABAAAB},\quad \overline{AAABAB},\quad \overline{AABABA},\quad \overline{BA},\quad \overline{BAAABA},\quad \overline{BABAAA}.
```

```@example notes
fig = Figure(size=(700, 700))
embeddedtreeplot!(juliaaxis(fig[1, 1], title="715/2047"), embedded(715//2047))
fig
```

*Hubbard tree for the internal address 1 -> 2 -> 6 -> 11, angle 715/2047, with each branch point labelled by its periodic kneading block.*

## Example: 1/9 -- tree, tripod images and preimage sketches

**Example:** ``\tfrac{1}{9}``.

The tree drawn at the top right has a single branch point of valence four. Along the horizontal arm, from left to right, sit the points ``2``, the branch point, ``0 = 6``, ``5``, ``4`` and the endpoint ``\beta``; the vertical arm carries ``1`` above the branch point and ``3`` below it, ending at ``-\beta``.

The left half of the sheet is a sequence of small working sketches in three colours, tracking a coloured tripod and its images. Two of them are annotated:

```math
\text{`double pre-image zone'} \qquad \text{`orbit structure and preimages'}
```

The last sketch in the sequence carries a spray of arrows pointing inward and outward from the branch point.

```@example notes
fig = Figure(size=(600, 600))
embeddedtreeplot!(juliaaxis(fig[1, 1], title="1/9"), embedded(1//9))
fig
```

*Tree for the angle 1/9: a valence-four branch point with the orbit points 1-6 and the endpoints ``\beta``, ``-\beta``.*

```@example notes
# The four arms at the branch point of the 1/9 tree, coloured, and their
# images under z ↦ z² + c: the first-return map rotates the arms.
HC = embedded(1//9)
c = HC.parameter
α = only(n for n in keys(HC.htree.adj) if length(HC.htree.adj[n]) == 4)
cross = arms(HC, α)
fig = Figure(size=(800, 250))
for k in 0:3
    ax = treeaxis(fig[1, k + 1], title=k == 0 ? "arms" : "f^$k")
    legsplot!(ax, [foldl((z, _) -> z .^ 2 .+ c, 1:k; init=arm) for arm in cross])
end
fig
```

*Column of small working sketches: a coloured tripod and its successive images, annotated ‘double pre-image zone’ and ‘orbit structure and preimages’.*

```@example notes
# The same arms pulled back by z ↦ ±√(z - c): one preimage cross at the
# branch point itself, the other at its negative.
fig = Figure(size=(500, 300))
ax = treeaxis(fig[1, 1])
legsplot!(ax, cross; linewidth=4)
legsplot!(ax, [pullback(arm, c) for arm in cross])
legsplot!(ax, [-pullback(arm, c) for arm in cross])
fig
```

*Four further cross-and-loop sketches in the lower right, repeating the sequence in monochrome and then in colour.*

## Segment-wise refinement in place, and plotting

Each vertex has one image so each line segment has one image.

For each segment:

- identify its image segment using kneading sequences.
- replace it by preimages of the image segment. (will involve `path_sqrt`

This can be done in-place, no need to copy the edges.

plotting…

## Orbit portraits and pinches

``O = \{z_1, \ldots, z_p\}`` a periodic orbit.

``A_i`` --- all the angles of rays which land at ``z_i``.

``P = \{A_1, \ldots, A_p\}`` orbit portrait.

Call ``A_i`` a “pinch”.

## Conjugate angles, four preimages, and choice of major leaf

Every periodic angle has a conjugate. Together they have 4 preimages corresponding to 2 choices of major leaf. There are only 2 choices because the major leaf needs to span both regions?

``\underline{\text{But}}`` if we take forward images of the minor leaf then the major leaf is unique.

One place we have 2 choices is picking the angle or its conjugate for the diameter dividing ``A`` and ``B``.

## Backwards construction from zero; next steps checklist

A different idea would be to go backwards, starting at zero (two neighbors) and using preimages.

First up: straight edge Hubbard tree.

- calculate the parameter with spider alg.
- calculate critical orbit with forward orbits
- calculate branch points with angles and external rays of Julia set from parameter.

## Parabolic parameters of ``\lambda(1+z/2)^2`` versus superattracting parameters of ``z^2+c``

If

```math
f:\; z \mapsto \lambda\left(1+\frac{z}{2}\right)^{2}
```

has a parabolic critical cycle at ``\lambda = \lambda_{*}``, then

```math
p:\; z \mapsto z^{2}+c
```

has a superattracting cycle at ``c = \lambda_{*}/2``.

```math
\begin{aligned}
  f'(z) &= 2\lambda\left(1+\frac{z}{2}\right)\cdot\frac{1}{2} \\
        &= \lambda\left(1+\frac{z}{2}\right) \\
  f'(0) &= \lambda
\end{aligned}
```

```math
\text{parabolic} \;\to\; |\lambda| = 1
```

```math
p'(z) = 2z \qquad p'(0) = 0
```

This is not the right way to think.

The superattracting critical orbit is ‘shadowed’ by a repelling orbit. These are the points where the external rays land.

## Four regions labelled by ``(B/A', 0/1)`` and their boundary colours

There are 4 regions

```math
\begin{aligned}
  B, 0  &\;-\; \text{red} \\
  A', 0 &\;-\; \text{blue} \\
  A', 1 &\;-\; \text{green} \\
  B', 1 &\;-\; \text{orange}
\end{aligned}
```

``(\beta, 0)`` marks the boundary between orange and red, so will be red orange.

``(0, -\beta)`` marks the boundary between blue and green, so will be blue green.

## Step 1: labelling the whole graph with 1 and 0, and labelling neighbours

Step 1: label the whole graph with 1 and 0. For the points on the boundary we have to label the neighbors with 1 and 0.

We will go ahead and do all the labelling on the neighbors, if a vertex is fully in a region then all its neighbors will have the same label since all the gaps will belong to the same region.

How to label neighbors: don't wanna do that.

label the whole tree: boundary is first. TODO: calculate the boudary and plot it. done.

## Three implementation options for computing vertex labels / external angles

Option 1: write a function that calculates 1 or 0 for a single vertex and call it repeatedly on each vertex in an orbit

Option 2: In one routine calculate all the external angles of all the tree's vertices at once.

Option 3?

Pick a chirality convention: each neighbor of a vertex will represent the gap counter clockwise of it.

## Shared iteration scheme for rays and spiders

In both cases:

- iterate over (target, source)
- pick the correct preimage of source

   - For rays --- by gluing distance
   - For spiders --- by region

- send this correct preimage to target

   - For rays --- append to end
   - For spiders --- add to kneading char

Question for spiders --- should all the feet go first?

Understand propositions 4.1 and 4.2 of spiders

## Four angles of an oriented Hubbard tree; spiders, legs and rays

Each oriented Hubbard tree has 4 angles

```math
\theta \qquad \bar{\theta} \qquad \theta^{*} \qquad \bar{\theta}^{*}
```

conjugate as well as swap 1s and zeros.

A spider is a sequence of legs.

A leg is a curve in the plane which goes to infinity and has a character to know which Domain it is in.

A ray is a curve in the plane which goes to infinity. (both are non-self intersecting)

## Node/vector signature reduction and ray counts by degree

```math
01011001011
```

```math
\texttt{function((node, vect\{node\}), (node, set\{node\}))}
```

```math
\Downarrow
```

```math
\texttt{(node, vect\{node\})}
```

A node with ``n`` neighbors could have

A node with degree ``n`` has how many external rays?

Critical point is a leaf and will have 2.

A branch point will have ``n``.

## Orienting the preimage of a branch point

Now we have an oriented branch point and we want to orient a preimage of this branch point.

target is oriented, we want to orient source.

- iterate over neighbors of source.

   - assign each neighbor the index of the neighbor of target towards the preimage of the neighbor

```@example notes
# In the 1/11 tree, a branch point off the periodic cycle (source) and the
# branch point it maps to (target).
HC = embedded(1//11)
adj = HC.htree.adj
image(n) = only(m for m in keys(adj) if m == Mandelbrot.shift(n))
isbranch(n) = length(adj[n]) >= 3
source = first(n for n in keys(adj) if isbranch(n) && isbranch(image(n)) &&
    !(n in orbit(image(n)).items))
target = image(source)
fig = Figure(size=(600, 600))
ax = juliaaxis(fig[1, 1], title="1/11")
embeddedtreeplot!(ax, HC; show_rays=false)
for (n, name) in ((source, "source"), (target, "target"))
    z = HC.vertices[n]
    scatter!(ax, [real(z)], [imag(z)], color=:purple, markersize=16)
    text!(ax, real(z), imag(z); text=name, offset=(8, -16), color=:purple)
end
fig
```

*Orienting a preimage: source branch point (right) maps to target branch point (left)*

## Orienting a Hubbard tree: cyclic order at characteristic branch points

Input: Hubbard tree

Output: Hubbard tree (oriented)

To embed the tree we will have to choose cyclic order for the characteristic branch points. The order of the remaining tree will follow.

Take a characteristic point of period ``k``. Consider the ``k``th return map.

Why must it permute the arms transitively?

Assuming it does permute transitively, the angle of the internal address tells us how much it rotates, or where the arm towards ``1`` is relative to the arm towards ``0``.

## Goal: interior and exterior coordinates; which disk coordinates to transport

Goal: calculate interior and exterior coordinates.

Prereq: what are the coordinates on the disk that we want to transport?

- angular coordinates: rays
- radial coordinates: equipotential lines

```math
\ln|z|
```

interior: pick ``\varepsilon`` and plot

```@example notes
fig, ax = expansionwedgeplot(3; size=(400, 400))
fig
```

*A disk cut by eight evenly spaced radii, drawn beside the note on angular and radial coordinates.*

## The data of the body; how to calculate a new body

The data of the body!

How to calculate a new Body:

1. find the two inverse images of ``\theta_1``. This will consist of a halving and a rotation depending on the angle of ``\theta_2``.
2. determine which sector is ``A`` by tracing the inverse image of leg ``2`` out to … the body. We know which one is correct because its foot will be at ``0``.

```@example notes
# The body for θ₁ = 9/31: its forward orbit θ₁, …, θ₅ on the circle, and the
# diameter through the preimages of θ₁ separating A from B.
θ = 9//31
fig = Figure(size=(450, 450))
partitiondisk!(treeaxis(fig[1, 1]), θ;
    marks=[mod(2^(k - 1) * θ, 1) => "θ$(Char(0x2080 + k))" for k in 1:5])
fig
```

*The body: a disk with five marked boundary angles and a chord separating regions A and B.*

```@example notes
fig = Figure(size=(250, 250))
partitiondisk!(treeaxis(fig[1, 1]), θ)
fig
```

*Small unlabelled disk with a single chord, drawn between steps 1 and 2.*

## Step 3: calculating the remaining legs by angle halving at the body

3. we can now calculate the remaining legs from the body down to the feet by doing angle halving at the body.

```@example notes
fig = Figure(size=(400, 400))
partitiondisk!(treeaxis(fig[1, 1]), θ; marks=[θ / 2 => L"\tilde{\theta}_1"])
fig
```

*Disk cut by a chord into sectors A and B, with the marked boundary angle tilde-theta-1.*

## Two new methods for spider legs; questions about Hubbard--Schleicher

So far, I've come up with two distinct methods of calculating spider legs distinct from the method provided by Hubbard and Schleicer.

H & S use intersections between the foot of the second leg (crossed out on the page)

H & S: what sector is

You need to understand the method they provide better. In particular the sentence

“we don't have these inverse images available”

and

“the lift of this segment starting at ``0`` leads to ``\tilde{w}``;”

## Angle for internal address [1,6,10]

What angle corresponds to ``[1,6,10]``?

```math
\text{Answer: } \tfrac{17}{1023} \approx \tfrac{4}{341}
```

Found by

```math
\begin{aligned}
A &\to 0\\
B &\to 1
\end{aligned}
```

```@example notes
# [1,6,10]: 1_{1/6} → 6_{1/2} → 10, which is 17/1023.
fig = Figure(size=(600, 600))
embeddedtreeplot!(juliaaxis(fig[1, 1], title="17/1023"), embedded(17//1023))
fig
```

*Tree of straight segments radiating from a small open circle, with numbered endpoints 0--9 and one edge inked green and labelled ``T_0``.*

## Denominators ``2^n \pm 1``: lollipop and star; embedding count is piecewise constant

```math
\begin{aligned}
\frac{1}{2^n+1} &\to \text{lollipop}\\[4pt]
\frac{1}{2^n-1} &\to \text{star}
\end{aligned}
```

The number of embeddings does not change as the denominator increases until we pass ``2^n``.

## Tree with second- and third-level branching (labels 0--3, 18, 19)

*(No prose on this page --- a single tree diagram only; see the figure.)*

```@example notes
# The page gives no angle. A six-armed branch point with critical-orbit
# labels up to 19 matches the period-20 angle 1/55 of the 1/6-limb.
fig = Figure(size=(700, 700))
embeddedtreeplot!(juliaaxis(fig[1, 1], title="1/55"), embedded(1//55); show_rays=false)
fig
```

*Tree of straight segments with three small open-circle vertices and filled-dot endpoints, labelled 0, 1, 2, 3, 18, 19.*

## Guessing the tree for ``1/17``

Make a guess for ``1/17``

``2``

A first attempt at the tree (see figure), with the four spokes labelled ``1/4``, ``2/3``, ``3/2``, ``4/1`` and the arm running out through the vertices ``0``, ``7``, ``6``, ``5``.

A second drawing of the same tree with the spokes in a different cyclic order: ``3/2``, ``1/4``, ``4/1``, ``2/3``, again with the arm ``0``, ``7``, ``6``, ``5``. Beside the arm: ``4?``

correct!

same for ``1/21``

```@example notes
# The 1/17 tree with no cyclic order chosen ...
fig = Figure(size=(700, 350))
treeplot!(treeaxis(fig[1, 1], title="1/17, unembedded"), 1//17; layout=:rt)
fig
```

*First guessed tree for ``1/17``: a five-valent hub with four short labelled spokes and one long arm*

```@example notes
# ... and as it actually embeds.
fig = Figure(size=(600, 600))
embeddedtreeplot!(juliaaxis(fig[1, 1], title="1/17"), 1//17)
fig
```

*Corrected tree for ``1/17``, same shape with the spokes in a different cyclic order*

## Angled internal address from an angle: allowed cyclic orders

How do we get an ``\underline{\text{angled}}`` internal address from an angle?

``1234``

``1324``

```math
\frac{4!}{4} = 3! = 6 \quad \text{orderings up to cyclic permutation}
```

Why are there only ``2`` allowable embeddings of this tree and not ``6``?

An ``\times``-shaped four-valent vertex with arms labelled ``3``, ``1``, ``2``, ``4`` (see figure) is marked *not allowed!*

map should be orientation preserving? cyclic order respecting

the two:

Two further ``\times``-shaped vertices, with arms labelled ``2``, ``1``, ``3``, ``4`` and ``4``, ``1``, ``3``, ``2`` (see figure).

```@example notes
# No parameter realises this cyclic order, so nothing in BrotViz can embed it;
# it is drawn by hand: arms 3, 1 at the top, 2, 4 at the bottom.
fig = Figure(size=(300, 300))
ax = treeaxis(fig[1, 1], title="not allowed")
for (φ, label) in ((3π/4, "3"), (π/4, "1"), (-3π/4, "2"), (-π/4, "4"))
    lines!(ax, [0, cos(φ)], [0, sin(φ)], color=:black)
    text!(ax, 1.15cos(φ), 1.15sin(φ); text=label, align=(:center, :center), fontsize=20)
end
scatter!(ax, [0], [0], color=:black)
fig
```

*The disallowed embedding: a four-valent vertex with arms labelled 3, 1, 2, 4*

```@example notes
# The two allowed embeddings of the four-armed star: the 1/4- and 3/4-limbs.
fig = Figure(size=(800, 400))
for (j, ang) in enumerate((1//4, 3//4))
    aia = AngledInternalAddress([1, 4], [ang])
    embeddedtreeplot!(juliaaxis(fig[1, j], title=addresstitle(aia)), aia)
end
fig
```

*The two allowed embeddings, drawn side by side*

## ``0.\overline{00001}``, and the kneading sequence of ``1/31``

What is ``0.\overline{00001}``?

```math
\frac{1}{32} \cdot \frac{1}{1 - \frac{1}{32}} = \frac{1}{31}
```

What is the kneading sequence of ``1/31``?

```math
a = 0.5 \qquad b = 16
```

```math
\begin{array}{ccccccccc}
1 & \to & 2 & \to & 4 & \to & 8 & \to & 16 \\
A & & A & & A & & A & & *
\end{array}
```

Something feels wrong. I don't think ``1/31`` & ``9/31`` should have the same kneading sequence.

*why not?*

*Is this what the angles in the internal address are about*

## Checking the doubling orbit of ``9/31``; generalizing

```math
a = 4.5 \qquad b = 20
```

check!

```math
\frac{9}{31} \to \frac{18}{31} \to \frac{5}{31} \to \frac{10}{31} \to \frac{20}{31}
```

```math
\begin{array}{ccccc}
0 & 1 & 0 & 0 & 1
\end{array}
```

(an arrow returns from ``20/31`` back to ``9/31``; in the margin, ``\frac{31}{2} = 15.5``)

*``\underline{\text{Generalize}}``*

*Given a binary expansion of period ``k``, calculate ``\dfrac{1}{1-2^{k}}`` then*

*``\underline{\text{TODO}}`` this but carefully and write it up*

What's the Hubbard tree for ``9/31``? First, what's its kneading sequence? ``AAAA*``!

```@example notes
# The doubling orbit of 9/31 on the circle. The upper half is digit 0 and the
# lower half digit 1, so the orbit reads off 0.(01001).
fig, ax = expansionwedgeplot(1; size=(450, 450))
for k in 0:4
    t = mod(2^k * 9//31, 1)
    scatter!(ax, [2π * t], [1.0], color=:black, markersize=12)
    text!(ax, 2π * t, 1.25; text="$(numerator(t))/31 → $(Int(t >= 1//2))",
        align=(:center, :center))
end
rlims!(ax, 0, 1.5)
fig
```

*The doubling orbit of ``9/31`` with its binary digits and a return arrow*

## Critical orbit for ``c=-1`` and embedding star graphs in the plane

Simple first…

Fix ``c = -1``, ``f(z) = z^2 - 1``:

```math
0 \longleftrightarrow -1
```

(everybody got ``\underline{2}`` preimages unless you are the critical value)

Our critical orbit

```math
0 \to c \to c^2 + c \to \cdots \to 0
\qquad\text{or}\qquad
0 \to c \to \cdots
```

has to be embedded in the plane in a way which is compatible with a continuous ``f``.

What do all the star graphs correspond to? All ‘A’ kneading sequences.

basilica ``\tfrac{1}{2}`` rabbit

```@example notes
# c = -1 is the centre of the 1/3 component: the basilica, 0 ↔ -1.
fig = Figure(size=(500, 400))
embeddedtreeplot!(juliaaxis(fig[1, 1], title="c = -1"), 1//3)
fig
```

*The two-cycle ``0 \leftrightarrow -1`` for ``c=-1``*

```@example notes
# A critical 3-cycle 0 → c → c² + c → 0, in its two embeddings.
fig = Figure(size=(800, 400))
for (j, ang) in enumerate((1//3, 2//3))
    aia = AngledInternalAddress([1, 3], [ang])
    embeddedtreeplot!(juliaaxis(fig[1, j], title=addresstitle(aia)), aia; show_rays=false)
end
fig
```

*The critical orbit drawn as a cycle, two versions*

```@example notes
fig = Figure(size=(800, 400))
embeddedtreeplot!(juliaaxis(fig[1, 1], title="basilica"), 1//3; show_rays=false)
embeddedtreeplot!(juliaaxis(fig[1, 2], title="rabbit"), 1//7; show_rays=false)
fig
```

*Basilica and rabbit star graphs with the basilica's angle 1/2*

## The rabbit's external angle: computing ``0.\overline{01001}``

The rabbit does not have external angle ``1/3``.

What angle is ``0.\overline{01001}``? ``0.0100101001``

```math
\begin{aligned}
&\tfrac{1}{4} + \tfrac{1}{32} \\
&\tfrac{1}{2^2} + \tfrac{1}{2^5} + \tfrac{1}{2^7} + \tfrac{1}{2^{10}} \\[2pt]
&\tfrac{1}{4}\left( \sum_{n=0}^{\infty} \left( \tfrac{1}{2^5} \right)^{n} \right) + \tfrac{1}{2^5}\left( \sum_{n=0}^{\infty} \left( \tfrac{1}{2^5} \right)^{n} \right) \\[2pt]
&\tfrac{1}{4}\left( \frac{1}{1 - \tfrac{1}{2^5}} \right) + \tfrac{1}{32}\left( \frac{1}{1 - \tfrac{1}{2^5}} \right) \\[2pt]
&\tfrac{1}{4}\cdot\tfrac{32}{31} + \tfrac{1}{32}\cdot\tfrac{32}{31} \;=\; \tfrac{9}{31}
\end{aligned}
```

```@example notes
fig = Figure(size=(400, 300))
labelledtreeplot!(treeaxis(fig[1, 1], title="AA*"), HubbardTree(KneadingSequence(1//7)); labels=:orbit)
fig
```

*Three-armed star graph with kneading sequence ``AA*``*

## Specification of `global_arms` and the local-arm sorting procedure

`Global arms(E, point)` will return a vector of vectors so that the first entry of each is the local arm and the rest follows the remaining points in the global arm.

---

For vertex in local arm:

- find the `preimage(global arms)`

Make table (preimage position, vertex).

Sort by preimage position.

Now cyclicly permute so the arm towards zero is last.

## Worked adjacency-list construction for a degree-4 vertex ``\alpha``

More concrete: say we have ``\alpha``.

1. Calculate the global arms of ``\alpha``.
2. ``\alpha`` is degree ``4`` so there are ``4`` global arms.
3. Create an entry for ``\alpha`` in the adjacency list

   ```math
   \alpha \quad [\; ,\; ,\; ,\; ]
   ```

Which of the global arms contains ``0``? The index of the “top” of that tree goes last:

```math
[\; ,\; ,\; ,\; a_0 ]
```

```math
[\, a_1,\, a_2,\, a_3,\, a_0 \,]
```

Now look at a preimage of ``\alpha`` and take the local arms. Look at the global arms of ``\alpha``.

## Arm positions at a characteristic point and inheriting order at preimages

Start -- characteristic point. The period of this point's orbit is ``k``.

```math
\begin{aligned}
&\text{the arm towards } 1 \text{ is in position } n && \text{numerator} = n\\
&\text{the arm towards } 1+k \text{ is in position } 2n\\
&\text{there are } d \text{ arms.} && \text{denominator} = d
\end{aligned}
```

The arm towards

```math
\begin{array}{cc}
1 & n\\
1+k & 2n\\
1+2k & 3n\\
\vdots & \vdots
\end{array}
```

Now take a preimage of this point, and the local arms of the preimage. For each arm of the preimage, find the global arm of the point that the local arm will be mapped to. Inherit the order from them.

## A shot at a tree implementation: Vector{Vector{Int}} versus a dictionary

A shot at a tree implementation

`Vector{Vector{Int}}`

```math
\begin{array}{l}
\text{index}\\
1\\
2\\
3\\
\vdots
\end{array}
```

Instead we want a dictionary which will allow us to work with subtrees.

Actually, do we ever really need the tree structure of global arms or do we just need to know the indices belonging to each? In the former case we can probably stick with the vector.

## Open questions: order of non-characteristic branch points, ``\alpha`` as root

How do we use dynamics to determine the order of branch points other than characteristic points?

Perhaps ``\alpha`` is good choice of root.

Question: is ``\alpha`` always in the set of marked points?

Question: is there a way to reconstruct the angle in the real unimodal case?

```math
c \in \mathbb{R}
```

## Rational angles, superattracting cycles and binary expansions

```math
z \mapsto z^2 + c
```

**Rational angle:**

```math
0 \le p/q < 1, \qquad p \text{ and } q \text{ coprime}
```

2 types:

```math
\left\{
\begin{array}{ll}
q \text{ Even} & - \text{ preperiodic} \\
q \text{ Odd} & - \text{ periodic}
\end{array}
\right.
```

Page 107, Corollary 10.4.10

Every quadratic polynomial with a superattracting cycle is conjugate to ``p_\theta`` for some ``\theta \in \mathbb{Q}^{\mathrm{odd}}/\mathbb{Z}``

``p_\theta`` is the quadratic polynomial equivalent to ``f_\theta``, the Thurston spider map defined by the ``\theta`` spider

binary ``\to`` angle

The binary expansion is a string of ones and zeros

```math
x_i, \quad i \in \mathbb{N}, \quad x_i \in \{0,1\}
```

```math
\theta = \sum_i x_i \cdot 2^{-i}
```

We describe rational angles with finite strings, consisting of a preperiod and a periodic part

```math
\theta = \sum_{i=1}^{\ell} x_i 2^{-i} + \sum_{i=\ell+1}^{\infty} x_i 2^{-i}, \qquad x_{i+k} = x_i \ \text{for } i > \ell
```

## Geometric series for the periodic tail; Hubbard tree questions

```math
\begin{aligned}
& x_{\ell+1} \cdot 2^{-(\ell+1)} \left( 1 + 2^{-k} + 2^{-2k} + \ldots \right) \\
+\ & x_{\ell+2} \cdot 2^{-(\ell+2)} \left( 1 + 2^{-k} + 2^{-2k} + \ldots \right)
\end{aligned}
```

```math
x_{\ell+1} \cdot 2^{-(\ell+1)} \cdot \frac{1}{1 - 2^{k}}
```

| Preperiod | Period |
|:--|:--|
| ``\displaystyle \sum_{i=1}^{\ell} x_i \cdot 2^{-i}`` | ``\displaystyle \sum_{i=\ell+1}^{\ell+k} x_i \cdot 2^{-i} \cdot \frac{1}{1 - 2^{-k}}`` |
|  | ``\displaystyle \frac{1}{1 - 2^{-k}} \sum_{i=\ell+1}^{\ell+k} x_i \cdot 2^{-i}`` |

---

**Hubbard trees**

What is a characteristic branch point? Given a Hubbard tree, how many characteristic branch points are there?

The characteristic branch point separates the critical value from zero

Given a Hubbard Tree, how many branch cycles are there? (last line runs off the bottom of the page)

