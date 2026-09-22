using BrotViz
using Test

using BrotViz: EscapeTimeProblem, escapetime, escape, triangulable, halton,
               banddisagreement, rawdisagreement
using Makie: Point2f
using BrotViz: smoothcount, value, build_mesh, planefit_prediction,
               interpolation_residuals, sample_arrays, DelTri, PENDING_VALUE

@testset "BrotViz.jl" begin

    @testset "samples" begin
        R2 = 100.0

        # The reference kernel. `escapetime` tests z *before* iterating and
        # reports the 1-based loop index, so with maxiter+1 it inspects
        # z_0 … z_maxiter and returns n+1 where `advance!` reports n. z_0 = 0
        # never escapes, so the extra check is free.
        function reference(c, maxiter)
            escapetime(EscapeTimeProblem(0.0im, z -> z * z + c, escape(R2), maxiter + 1))
        end

        cases = [0.0 + 0.0im, -1.0 + 0.0im, -0.123 + 0.745im, 0.25 + 0.0im,
                 1.0 + 0.0im, 0.5 + 0.0im, -2.1 + 0.0im, 0.3 + 0.5im,
                 -0.75 + 0.1im, -1.75 + 0.02im, 2.0 + 2.0im]

        @testset "advance! agrees with escapetime (chunk = $chunk)" for chunk in (1, 7, 1000)
            for maxiter in (64, 512)
                for c in cases
                    s = BrotViz.Sample(c)
                    while !BrotViz.advance!(s, chunk, R2, maxiter)
                    end
                    iter, z = reference(c, maxiter)
                    if s.status === Escaped
                        @test iter == s.n + 1
                        @test z ≈ s.z
                        @test s.n <= maxiter
                    else
                        @test s.status === Interior
                        @test isnan(iter)
                        @test s.n == maxiter
                    end
                end
            end
        end

        @testset "chunk size does not change the answer" begin
            for c in cases
                results = map((1, 7, 1000)) do chunk
                    s = BrotViz.Sample(c)
                    while !BrotViz.advance!(s, chunk, R2, 512)
                    end
                    (s.n, s.status)
                end
                @test allequal(results)
            end
        end

        @testset "round-robin fairness" begin
            set = SampleSet(; escape_radius2 = R2, maxiter = 10_000)
            for c in cases
                push_sample!(set, c)
            end
            chunk = 5
            for k in 1:40
                run_round!(set, chunk)
                # every sample still pending has had exactly the same work
                for idx in set.pending
                    @test set.samples[idx].n == k * chunk
                end
            end
            @test npending(set) > 0   # the interior points are still going
        end

        @testset "pending bookkeeping" begin
            set = SampleSet(; escape_radius2 = R2, maxiter = 200)
            for c in cases
                push_sample!(set, c)
            end
            @test npending(set) == length(cases)
            @test nresolved(set) == 0
            for chunk in (1, 1, 3, 17, 200)
                run_round!(set, chunk)
                expected = Set(i for i in eachindex(set.samples)
                               if set.samples[i].status === Pending)
                @test Set(set.pending) == expected
                @test length(set.pending) == length(expected)   # no duplicates
                @test nresolved(set) == length(set) - length(expected)
            end
            @test isempty(set.pending)
            @test all(s -> s.status !== Pending, set.samples)
        end

        @testset "known classifications" begin
            set = SampleSet(; escape_radius2 = R2, maxiter = 512)
            interior = [0.0 + 0.0im, -1.0 + 0.0im, -0.123 + 0.745im]
            exterior = [1.0 + 0.0im, 0.5 + 0.0im, -2.1 + 0.0im]
            for c in vcat(interior, exterior)
                push_sample!(set, c)
            end
            resolve_all!(set)
            for (i, _) in enumerate(interior)
                @test set.samples[i].status === Interior
                @test set.samples[i].n == 512
            end
            for (j, _) in enumerate(exterior)
                s = set.samples[length(interior) + j]
                @test s.status === Escaped
                @test 0 < s.n < 20     # these leave the disc of radius 10 fast
            end
        end

        @testset "dedupe on insert" begin
            set = SampleSet()
            @test push_sample!(set, 0.1 + 0.1im; mindist = 0.05) == 1
            @test push_sample!(set, 0.11 + 0.1im; mindist = 0.05) === nothing
            @test push_sample!(set, 0.3 + 0.1im; mindist = 0.05) == 2
            @test length(set) == 2
            @test push_sample!(set, 0.11 + 0.1im) == 3   # no tolerance, always adds
        end

        @testset "empty!" begin
            set = SampleSet()
            sprinkle!(set, (-2, 1), (-1, 1), 20)
            @test length(set) == 20
            empty!(set)
            @test isempty(set)
            @test npending(set) == 0
        end
    end

    @testset "samplers" begin
        @test halton(1, 2) ≈ 1 / 2
        @test halton(2, 2) ≈ 1 / 4
        @test halton(3, 2) ≈ 3 / 4
        @test halton(1, 3) ≈ 1 / 3
        @test halton(2, 3) ≈ 2 / 3

        set = SampleSet()
        @test sprinkle!(set, (-2.0, 1.0), (-1.0, 1.0), 200) == 200
        @test all(s -> -2.0 <= real(s.c) <= 1.0 && -1.0 <= imag(s.c) <= 1.0, set.samples)
        # a second sprinkle continues the sequence rather than repeating it
        sprinkle!(set, (-2.0, 1.0), (-1.0, 1.0), 200)
        @test length(unique(s.c for s in set.samples)) == 400
    end

    @testset "refine!" begin
        XL, YL = (-2.5, 1.0), (-1.25, 1.25)

        fresh(n) = begin
            set = SampleSet(; maxiter = 256)
            sprinkle!(set, XL, YL, n)
            resolve_all!(set)
            set
        end

        @testset "degenerate inputs are no-ops" begin
            @test refine!(SampleSet(), 10) == (added = 0, removed = 0)   # empty
            set = SampleSet(; maxiter = 64)
            push_sample!(set, 0.0 + 0.0im); resolve_all!(set)
            @test refine!(set, 10) == (added = 0, removed = 0)   # one resolved sample
            # three collinear samples cannot be triangulated
            set = SampleSet(; maxiter = 64)
            for c in (-1.0 + 0im, 0.0 + 0im, 0.5 + 0im)
                push_sample!(set, c)
            end
            resolve_all!(set)
            @test refine!(set, 10) == (added = 0, removed = 0)
            # pending samples have no value yet, so they do not count
            set = fresh(40)
            sprinkle!(set, XL, YL, 40; offset = 1000)       # 40 fresh pending
            before = length(set)
            @test refine!(set, 5; prune = false).added == 5
            @test length(set) == before + 5
        end

        @testset "unanimous samples have nothing to refine" begin
            # every sample deep inside the set resolves to Interior at maxiter,
            # so every edge has zero disagreement and the walk stops at once
            set = SampleSet(; maxiter = 64)
            sprinkle!(set, (-0.3, -0.1), (-0.1, 0.1), 50)
            resolve_all!(set)
            @test all(s -> s.status === Interior, set.samples)
            @test refine!(set, 20; prune = false) == (added = 0, removed = 0)
        end

        @testset "adds pending midpoints inside the existing hull" begin
            set = fresh(200)
            xs = [real(s.c) for s in set.samples]
            ys = [imag(s.c) for s in set.samples]
            before = length(set)
            added = refine!(set, 50; prune = false).added
            @test added == 50
            @test length(set) == before + 50
            new = set.samples[(before + 1):end]
            @test all(s -> s.status === Pending, new)
            @test all(s -> s.n == 0, new)
            # a midpoint of two existing samples cannot escape their bounds
            @test all(s -> minimum(xs) <= real(s.c) <= maximum(xs), new)
            @test all(s -> minimum(ys) <= imag(s.c) <= maximum(ys), new)
        end

        @testset "mindist rejections do not consume budget" begin
            set = fresh(200)
            before = length(set)
            # a tolerance far larger than the sample spacing rejects everything
            @test refine!(set, 30; mindist = 10.0, prune = false).added == 0
            @test length(set) == before
            @test refine!(set, 30; mindist = 1.0e-9, prune = false).added == 30
        end

        @testset "within restricts refinement to a window" begin
            set = fresh(3000)
            win = ((-1.2, -0.4), (0.0, 0.6))     # straddles ∂M
            before = length(set)
            added = refine!(set, 60; mindist = 1.0e-6, within = win, prune = false).added
            @test added == 60
            new = set.samples[(before + 1):end]
            # midpoints of edges between in-window samples stay in the window
            @test all(new) do s
                win[1][1] <= real(s.c) <= win[1][2] && win[2][1] <= imag(s.c) <= win[2][2]
            end
            # a window holding fewer than three samples is a no-op
            @test refine!(set, 10; within = ((0.99, 1.0), (1.2, 1.25))) ==
                  (added = 0, removed = 0)
        end

        @testset "budget is capped by the number of disagreeing edges" begin
            # A window lying almost entirely inside the set has only a handful
            # of edges with anything to say, and refine! stops there rather
            # than inserting dead points to fill the quota. Asking again after
            # the new samples resolve is what makes progress — which is
            # exactly what pressing `r` repeatedly does.
            set = fresh(400)
            win = ((-1.0, -0.5), (0.0, 0.4))
            first_pass = refine!(set, 1000; criterion = :disagreement,
                                 mindist = 1.0e-6, within = win, prune = false).added
            @test 0 < first_pass < 1000
            resolve_all!(set)
            @test refine!(set, 1000; criterion = :disagreement,
                          mindist = 1.0e-6, within = win, prune = false).added > 0
        end

        @testset "priority policies rank differently" begin
            # raw disagreement ignores length, band disagreement does not, so
            # they must not pick the same midpoints
            a = fresh(200); b = fresh(200)
            refine!(a, 40; priority = rawdisagreement, prune = false)
            refine!(b, 40; priority = banddisagreement, prune = false)
            pa = Set(s.c for s in a.samples[201:end])
            pb = Set(s.c for s in b.samples[201:end])
            @test pa != pb
        end

        @testset "pruning leaves the picture identical" begin
            # The whole justification for deleting samples: a sample whose
            # neighbours all agree is invisible, so the nearest-neighbour
            # reconstruction must come out bit-identical afterwards.
            set = fresh(1500)
            G = 140
            # offset off the sample lattice so exact equidistant ties, which
            # could legitimately flip, do not arise
            gx = range(XL[1] + 0.0037, XL[2] - 0.0037; length = G)
            gy = range(YL[1] + 0.0029, YL[2] - 0.0029; length = G)
            function reconstruct(set)
                sx = [real(s.c) for s in set.samples]
                sy = [imag(s.c) for s in set.samples]
                sv = [s.status === Interior ? -1 : s.n for s in set.samples]
                [begin
                     best, bd = 1, Inf
                     for k in eachindex(sx)
                         d = (sx[k] - x)^2 + (sy[k] - y)^2
                         d < bd && (bd = d; best = k)
                     end
                     sv[best]
                 end for x in gx, y in gy]
            end

            before_img = reconstruct(set)
            before_n = length(set)
            r = refine!(set, 0; mindist = 1.0e-6)     # prune only, add nothing
            @test r.added == 0
            @test r.removed > 0
            @test length(set) == before_n - r.removed
            @test reconstruct(set) == before_img
        end

        @testset "pruning strips a flat region to its rim" begin
            # every sample deep inside the set resolves Interior, so only the
            # convex hull — the samples whose cells reach the clip rect — has
            # any reason to survive
            set = SampleSet(; maxiter = 64)
            sprinkle!(set, (-0.3, -0.1), (-0.1, 0.1), 200)
            resolve_all!(set)
            @test all(s -> s.status === Interior, set.samples)
            r = refine!(set, 100)
            @test r.added == 0
            @test length(set) < 20            # just the hull
            @test length(set) >= 3            # still triangulable
            # and there is nothing left to do
            @test refine!(set, 100) == (added = 0, removed = 0)
        end

        @testset "pruning keeps pending bookkeeping straight" begin
            # deletion shifts every index past it, so `pending` has to be
            # rebuilt — this is the test that catches getting that wrong
            set = fresh(600)
            sprinkle!(set, (-0.6, 0.3), (-0.5, 0.5), 150; offset = 9000)
            @test npending(set) > 0
            pending_before = Set(set.samples[i].c for i in set.pending)

            r = refine!(set, 20; mindist = 1.0e-6)
            @test r.removed > 0

            expected = Set(i for i in eachindex(set.samples)
                           if set.samples[i].status === Pending)
            @test Set(set.pending) == expected
            @test length(set.pending) == length(expected)
            # pending samples are never deleted (they have no value to agree on)
            @test pending_before ⊆ Set(set.samples[i].c for i in set.pending)
            # and the set still runs correctly afterwards
            resolve_all!(set)
            @test isempty(set.pending)
            @test all(s -> s.status !== Pending, set.samples)
        end

        @testset "prunetol widens what counts as agreement" begin
            set = fresh(1500)
            loose = refine!(deepcopy(set), 0; prunetol = 5.0).removed
            exact = refine!(deepcopy(set), 0; prunetol = 0.0).removed
            @test loose > exact > 0
        end

        @testset "refinement beats a plain sprinkle at the boundary" begin
            # The whole claim of §2a.4: samples should concentrate on ∂M
            # without being told it exists. Measured against ground truth on a
            # grid, as misclassified area under nearest-neighbour
            # reconstruction. Fully deterministic — Halton and the edge
            # ranking are both reproducible.
            G = 100
            gx = range(XL[1], XL[2]; length = G)
            gy = range(YL[1], YL[2]; length = G)
            function isinterior(c)
                s = BrotViz.Sample(c)
                while !BrotViz.advance!(s, 256, 100.0, 256)
                end
                s.status === Interior
            end
            truth = [isinterior(complex(x, y)) for x in gx, y in gy]

            function misclassified(set)
                sx = [real(s.c) for s in set.samples]
                sy = [imag(s.c) for s in set.samples]
                si = [s.status === Interior for s in set.samples]
                bad = 0
                for (i, x) in enumerate(gx), (j, y) in enumerate(gy)
                    best, bd = 1, Inf
                    for k in eachindex(sx)
                        d = (sx[k] - x)^2 + (sy[k] - y)^2
                        d < bd && (bd = d; best = k)
                    end
                    si[best] == truth[i, j] || (bad += 1)
                end
                bad / length(truth)
            end

            budget = 800
            plain = fresh(budget)

            refined = fresh(200)
            while length(refined) < budget
                refine!(refined, min(200, budget - length(refined));
                        mindist = 1.0e-6, prune = false)
                resolve_all!(refined)
            end

            @test length(refined) == length(plain) == budget
            @test all(s -> s.status !== Pending, refined.samples)
            @test misclassified(refined) < misclassified(plain)
        end
    end

    @testset "smoothcount" begin
        R2 = 100.0
        resolve(c, maxiter = 512) = begin
            s = BrotViz.Sample(c)
            while !BrotViz.advance!(s, 64, R2, maxiter)
            end
            s
        end

        @test isnan(smoothcount(resolve(-0.2 + 0.0im), 64), )  skip = false
        @test isnan(smoothcount(resolve(-0.2 + 0.0im), 64))      # interior

        @testset "brackets its integer count" begin
            for c in (1.0 + 0im, 0.5 + 0.2im, -2.1 + 0im, 0.3 + 0.5im, -0.75 + 0.3im)
                s = resolve(c)
                ν = smoothcount(s, R2)
                @test s.n <= ν < s.n + 1
            end
        end

        @testset "is continuous where the integer count jumps" begin
            # Walk a ray outwards across a band boundary. The integer count
            # steps by one; ν must not.
            prev = nothing
            jumped = false
            maxjump = 0.0
            for t in range(0.60, 0.90; length = 4000)
                s = resolve(complex(t, 0.35))
                s.status === Escaped || continue
                ν = smoothcount(s, R2)
                if prev !== nothing
                    prev[1] == s.n || (jumped = true)
                    maxjump = max(maxjump, abs(ν - prev[2]))
                end
                prev = (s.n, ν)
            end
            @test jumped                # the integer count really does step
            @test maxjump < 0.05        # ν does not
        end

        @test value(SampleSet(), resolve(1.0 + 0im); smooth = false) == resolve(1.0 + 0im).n
    end

    @testset "sample_arrays" begin
        set = SampleSet(; maxiter = 256)
        for c in (-0.2 + 0.0im, 1.0 + 0.0im, 0.4 + 0.4im)
            push_sample!(set, c)
        end
        push_sample!(set, -1.0 + 0.6im)          # left pending
        resolve_all!(SampleSet())                 # no-op, keeps the last pending
        for i in 1:3
            while !BrotViz.advance!(set.samples[i], 256, set.escape_radius2, set.maxiter)
            end
        end
        filter!(!=(4), set.pending)

        pos, val, cr = sample_arrays(set; admit = :immediate, smooth = true)
        @test length(pos) == length(val) == 4
        @test val[4] == PENDING_VALUE             # pending routes to lowclip
        @test isnan(val[1])                       # interior routes to nan_color
        @test cr[1] >= 1                          # log10 colorscale stays in domain

        pos2, val2, _ = sample_arrays(set; admit = :onresolve, smooth = true)
        @test length(pos2) == length(val2) == 3   # pending excluded
        @test all(!=(PENDING_VALUE), val2)
    end

    @testset "build_mesh" begin
        pts = Point2f[(0, 0), (1, 0), (1, 1), (0, 1), (0.5, 0.4)]
        tri = DelTri.triangulate([(Float64(p[1]), Float64(p[2])) for p in pts];
                                 randomise = false)
        ntri = length(collect(DelTri.each_solid_triangle(tri)))

        @testset "all interpolatable: one face per triangle, no new vertices" begin
            vals = Float32[1, 2, 3, 4, 5]
            verts, faces, vcols = build_mesh(pts, vals, tri)
            @test verts == pts
            @test vcols == vals
            @test length(faces) == ntri
        end

        @testset "a discontinuity splits its triangles into flat quads" begin
            vals = Float32[1, 2, 3, 4, NaN32]      # the middle vertex is interior
            verts, faces, vcols = build_mesh(pts, vals, tri)
            @test length(verts) == length(vcols)
            @test length(verts) > length(pts)      # midpoints and centroids added
            @test verts[1:length(pts)] == pts      # originals keep their indices
            @test all(f -> all(i -> 1 <= i <= length(verts), f), faces)
            # every added vertex copies the colour of the sample that owns it,
            # so no face straddling the discontinuity ever interpolates
            straddling = count(t -> any(isnan, (vals[t[1]], vals[t[2]], vals[t[3]])),
                               DelTri.each_solid_triangle(tri))
            @test straddling > 0
            @test length(faces) == (ntri - straddling) + 6 * straddling
            @test length(verts) == length(pts) + 9 * straddling
        end

        @testset "pending vertices are treated as a discontinuity too" begin
            # you cannot Gouraud-interpolate into a value that does not exist yet
            a = build_mesh(pts, Float32[1, 2, 3, 4, PENDING_VALUE], tri)
            b = build_mesh(pts, Float32[1, 2, 3, 4, NaN32], tri)
            @test length(a[2]) == length(b[2])
        end
    end

    @testset "interpolation residuals" begin
        @testset "a plane is predicted exactly" begin
            # w = 3 + 2x - y sampled on an irregular neighbourhood
            dx = [1.0, -0.7, 0.2, -1.3, 0.9]
            dy = [0.3, 1.1, -0.9, -0.2, 0.7]
            w = [3 + 2dx[i] - dy[i] for i in eachindex(dx)]
            @test planefit_prediction(dx, dy, w) ≈ 3.0 atol = 1e-10
        end
        @testset "collinear neighbours fall back to the mean" begin
            dx = [1.0, 2.0, -1.0]; dy = [1.0, 2.0, -1.0]   # all on y = x
            w = [2.0, 4.0, 6.0]
            @test planefit_prediction(dx, dy, w) ≈ 4.0
        end
        @testset "classes" begin
            pts = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0), (0.4, 0.4)]
            nbrs = [[2, 3, 5], [1, 4, 5], [1, 4, 5], [2, 3, 5], [1, 2, 3, 4]]
            isint = [false, false, false, false, false]
            w = [0.0, 1.0, 1.0, 2.0, 0.8]            # exactly planar: w = x + y
            resid, span = interpolation_residuals(pts, w, isint, nbrs)
            @test resid[5] ≈ 0.0 atol = 1e-9         # 0.4+0.4 = 0.8, predicted
            @test span ≈ 2.0
            # interior among interior is free; a mixed neighbourhood is not
            isint2 = [true, true, true, true, true]
            r2, _ = interpolation_residuals(pts, w, isint2, nbrs)
            @test all(iszero, r2)
            isint3 = [false, false, false, false, true]
            r3, _ = interpolation_residuals(pts, w, isint3, nbrs)
            @test r3[5] > 0 && all(i -> r3[i] > 0, 1:4)
        end
    end

    @testset "triangulable" begin
        @test !triangulable(Point2f[])
        @test !triangulable(Point2f[(0, 0), (1, 0)])
        @test !triangulable(Point2f[(0, 0), (0, 0), (0, 0)])
        @test !triangulable(Point2f[(0, 0), (1, 1), (2, 2)])          # collinear
        @test !triangulable(Point2f[(0, 0), (1, 1), (2, 2), (3, 3)])
        @test !triangulable(Point2f[(0, 0), (0, 1), (0, 2)])          # vertical
        @test triangulable(Point2f[(0, 0), (1, 0), (0, 1)])
        @test triangulable(Point2f[(0, 0), (1, 1), (2, 2), (0, 2)])
        # scale invariance: the same shape, shrunk, is still triangulable
        @test triangulable(Point2f[(0, 0), (1e-4, 0), (0, 1e-4)])
    end
end
