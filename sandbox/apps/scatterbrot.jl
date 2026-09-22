# scatterbrot — click-to-sample Mandelbrot explorer (Phase 1 MVP)
#
# Pan, zoom, and click to place sample points in the parameter plane. A click
# queues a batch of samples, which tile the plane grey immediately and then
# colour in one cell at a time as the scheduler resolves them. Points advance
# in strict round-robin, so all of them get the same number of iterations per
# unit wall clock. The UI never blocks: the scheduler runs on the render
# thread in slices of a few milliseconds, so interaction latency is bounded
# by the budget rather than by the workload.
#
# Because the samples *are* complex numbers and the cells *are* polygons in
# data space, pan and zoom are exact at every scale and old data never goes
# stale. Zooming in does not invalidate anything — it just reveals that your
# cells are too coarse, which is the cue to add more points.
#
# ---------------------------------------------------------------------------
# Running it
# ---------------------------------------------------------------------------
# This script lives outside the package so GLMakie never enters
# BrotViz/Project.toml. Give it an environment of its own:
#
#   julia -e 'using Pkg; Pkg.activate("/tmp/scatterbrot"); Pkg.add("GLMakie");
#             Pkg.develop(path="'$PWD'")'
#   julia --project=/tmp/scatterbrot sandbox/apps/scatterbrot.jl
#
# Phase 1 is single threaded, so plain `julia` is enough. Phase 4 (threaded
# rounds) will want `julia -t auto,1` — this machine has 16 CPU threads but
# Julia still defaults to 1, and the `,1` reserves an interactive thread.
#
# ---------------------------------------------------------------------------
# Controls
# ---------------------------------------------------------------------------
#   left click      sprinkle a cluster of CLICK_N samples around the pointer
#   left drag       rectangle zoom (Makie's :rectanglezoom, untouched)
#   right drag      pan            (Makie's :dragpan, untouched)
#   scroll          zoom in / out  (Makie's :scrollzoom, untouched)
#   0  /  ctrl+click   back to the home view
#   s               sprinkle SPRINKLE_N Halton points over the whole viewport
#   r               refine once: add where the data disagrees, drop what is
#                   redundant
#   f               toggle continuous refinement at REFINE_HZ
#   c               clear every sample
#   t               toggle reconstruction: linear mesh / Voronoi cells
#   a               toggle cell admission: admit-immediately / paint-on-resolve
#   space           pause / resume the scheduler
#
# `0` and ctrl+click do the same thing. ctrl+click is Makie's own
# `:limitreset` interaction and comes for free — which is exactly why
# `:addsample` has to ignore clicks with ctrl held, or resetting the view
# would also dump a cluster of samples wherever you happened to click.
#
# Left click is safe to bind even though left drag is rectangle zoom: Makie's
# mouse state machine only emits `leftclick` on a press+release with no
# intervening drag, and emits leftdragstart/leftdrag/leftdragstop instead when
# there is one. So clicking and rectangle-zooming do not fight.
#
# ---------------------------------------------------------------------------
# Manual checklist (this script is not unit tested; src/samples.jl is)
# ---------------------------------------------------------------------------
#   [ ] window opens showing only axes, says it is warming up for ~15 s, then
#       says ready — and the first click after that lands in well under a second
#   [ ] with WARMUP = false the first click instead hangs for ~15 s (this is
#       what the warm-up is buying)
#   [ ] the first click alone brings up a tessellation (a cluster is enough),
#       grey at first, colouring in within a beat
#   [ ] a click at 15k+ samples repaints promptly — one retessellation for the
#       batch, not one per resolving sample
#   [ ] clicking rapidly places a cluster per click — none are swallowed as
#       double clicks
#   [ ] ctrl+click returns to the home view and does NOT add samples
#   [ ] `0` returns to the home view from any zoom
#   [ ] third non-collinear click brings up the tessellation
#   [ ] three clicks along the real axis do NOT crash; the fourth off-axis
#       click brings the tessellation up
#   [ ] clicking inside the main cardioid leaves a grey cell for a visible
#       beat, then resolves to black
#   [ ] right-dragging pans without placing a sample; cells move with the data
#   [ ] `s` fills the viewport with grey and the cells colour in, with the
#       last grey to survive tracing ∂M
#   [ ] `s` then `r` held down: samples visibly pile onto ∂M and onto the
#       escape-time band boundaries, and the open exterior stays sparse
#   [ ] `r` with nothing on screen is a no-op, not an error
#   [ ] `r` over a view dominated by the interior *drops* samples — the count
#       can go down, and the picture does not change when it does
#   [ ] `f` leaves it refining by itself; the status line shows auto-refine
#       and the sample count climbs and settles
#   [ ] `f` refinement keeps the window responsive to pan and zoom, and
#       visibly slows down as the sample count climbs
#   [ ] `space` freezes the remaining grey cells; pressing it again resumes
#   [ ] `c` empties the window without error, and clicking afterwards rebuilds
#   [ ] `a` flips between grey-pending cells and not drawing them at all
#       (and `:onresolve` is visibly the slower of the two at size)
#   [ ] `t` flips between the smooth mesh and the Voronoi cells, keeping the
#       samples; the mesh is visibly smoother in the exterior at equal count
#   [ ] the mesh leaves the corners of the view black until samples reach
#       them — it only covers the convex hull. `s` fixes it, clicks do not
#   [ ] zooming in while thousands of points are pending stays responsive

using GLMakie
using BrotViz
using BrotViz: limits_of

# --- knobs ------------------------------------------------------------------

const MAXITER      = 512            # low on purpose: without period detection
                                    # every interior point runs the full budget
                                    # and takes a slice of every round (§7.7)
const CHUNK        = 64             # iterations per sample per round
const IDLE_BUDGET  = 8_000_000      # ns per frame when the user is not touching it
const BUSY_BUDGET  = 2_000_000      # ns per frame while panning/zooming/clicking
const SPRINKLE_N   = 400
const REFINE_N     = 400
const REFINE_PAD   = 0.1            # grow the refine window by this much per side
const REFINE_HZ    = 4              # continuous-refinement rate (`f`), when
                                    # passes are cheap enough to keep it up
const REFINE_DUTY  = 0.25           # ...and never let refinement eat more than
                                    # this share of wall clock, whatever it costs
const REFINE_WAIT  = 2.0            # ...but never wait longer than this, so one
                                    # freak pass cannot stall the loop
const MAX_SAMPLES  = 30_000         # continuous refinement stops growing here;
                                    # past this the retriangulation per pass is
                                    # a visible hitch
const CLICK_N      = 60             # samples per click
const CLICK_FRAC   = 0.02           # cluster radius, as a fraction of the viewport
const RECON        = :mesh          # :mesh (piecewise linear, smooth exterior)
                                    # or :voronoi (piecewise constant, honest
                                    # about knowing nothing between samples)
const WARMUP       = true           # compile the draw path at startup rather
                                    # than during the user's first click
const COMMIT_HZ    = 30             # cap on republishing to the plot. Resolved
                                    # samples land within 33 ms, which is below
                                    # notice, and it halves the per-frame cost
                                    # of pushing colours at large sample counts
const DEDUPE_FRAC  = 0.002          # min sample spacing, as a fraction of the
                                    # viewport width (coincident generators make
                                    # DelaunayTriangulation warn and drop points)

# --- scene ------------------------------------------------------------------

set = SampleSet(; maxiter = MAXITER)

const HOME = (-2.5, 1.0, -1.25, 1.25)

fig = Figure(size = (1100, 900))
ax = Axis(fig[1, 1]; aspect = DataAspect(),
          xlabel = "Re c", ylabel = "Im c",
          title = "scatterbrot — click to sample")
limits!(ax, HOME...)

make_plot(kind) = kind === :mesh ? meshbrotplot!(ax, set) : voronoibrotplot!(ax, set)

recon = Ref(RECON)
vb = Ref{Any}(make_plot(RECON))

status = Observable("0 samples")
Label(fig[2, 1], status; tellwidth = false, halign = :left, font = :regular)

pad((lo, hi), frac) = (lo - frac * (hi - lo), hi + frac * (hi - lo))

# --- interaction ------------------------------------------------------------

paused = Ref(false)
auto = Ref(false)          # continuous refinement (`f`)
next_commit = Ref(0.0)
next_refine = Ref(0.0)     # earliest wall-clock time for the next auto pass
dirty = Ref(true)          # something changed that commit! has not published yet

# One refinement pass over what is on screen. Refining samples nobody can see
# would cost a bigger retriangulation for a result nobody is looking at, and
# the retriangulation is the whole cost.
function swap_reconstruction!()
    old = vb[]
    detach!(old)                     # drops the plot, keeps the samples
    recon[] = recon[] === :mesh ? :voronoi : :mesh
    vb[] = make_plot(recon[])
    vb[].admit = old.admit
    return nothing
end

function refine_view!()
    xlims, ylims = limits_of(ax)
    r = refine!(set, REFINE_N;
                mindist = DEDUPE_FRAC * (xlims[2] - xlims[1]),
                within = (pad(xlims, REFINE_PAD), pad(ylims, REFINE_PAD)))
    dirty[] = true
    return r
end

# `event.data` is already in data coordinates, so there is nothing to unproject.
#
# Both click types are handled on purpose. Makie's state machine emits
# `leftdoubleclick` *instead of* `leftclick` when a second click lands within
# 200 ms of the first, anywhere in the axis — so listening only for
# `leftclick` silently drops every other sample when you click quickly, which
# is exactly how you use this app.
register_interaction!(ax, :addsample) do event::MouseEvent, ax
    if event.type === MouseEventTypes.leftclick ||
            event.type === MouseEventTypes.leftdoubleclick
        # ctrl+click belongs to Makie's :limitreset
        ispressed(ax.scene, Keyboard.left_control) && return
        (xmin, xmax), _ = limits_of(ax)
        width = xmax - xmin
        c = complex(event.data[1], event.data[2])
        mindist = DEDUPE_FRAC * width
        # the sample exactly where they clicked, then a cluster around it —
        # one point is rarely what you meant, and a cluster gives refinement
        # something to chew on straight away
        push_sample!(set, c; mindist)
        sprinkle_disc!(set, c, CLICK_FRAC * width, CLICK_N - 1; mindist)
        dirty[] = true
    end
end

# Idleness signal for the frame budget: spend little while the user is
# actively manipulating the view, spend more when they are not. `tick.state`
# is useless for this — we dirty the scene ourselves every frame we do work,
# so a SkippedRenderTick never arrives while there is anything to compute.
#
# Both halves of this have to route around Makie's event *consumption*, and
# the obvious spellings are silently broken:
#
#   * Tracking the button with `on(events.mousebutton)` and a `pressed` flag
#     latches on forever the first time an interaction consumes a click.
#     ctrl+click does exactly that — it is Makie's `:limitreset` — so one
#     view reset would pin the budget at BUSY forever and stop auto-refine
#     dead, because the flag sees the press and never the release. Read
#     `mousebuttonstate` instead: it is the actual button state, maintained
#     upstream of any listener, and cannot go stale.
#   * `on(events.scroll)` at the default priority never fires at all, because
#     `:scrollzoom` consumes the event first. Ask for a higher priority and
#     return nothing, so we observe it without taking it away.
last_scroll = Ref(0.0)
on(events(ax.scene).scroll, priority = 100) do _
    last_scroll[] = time()
    return Consume(false)
end
manipulating() = !isempty(events(ax.scene).mousebuttonstate) ||
                 time() - last_scroll[] < 0.2

on(events(fig).keyboardbutton) do e
    e.action === Keyboard.press || return Consume(false)
    if e.key === Keyboard.c
        empty!(vb[])
        dirty[] = true
    elseif e.key === Keyboard.s
        xlims, ylims = limits_of(ax)
        sprinkle!(set, xlims, ylims, SPRINKLE_N;
                  mindist = DEDUPE_FRAC * (xlims[2] - xlims[1]))
        dirty[] = true
    elseif e.key === Keyboard.r
        refine_view!()
    elseif e.key === Keyboard.f
        auto[] = !auto[]
        dirty[] = true
    elseif e.key === Keyboard._0
        limits!(ax, HOME...)
        dirty[] = true
    elseif e.key === Keyboard.t
        swap_reconstruction!()
        dirty[] = true
    elseif e.key === Keyboard.a
        vb[].admit = vb[].admit === :onresolve ? :immediate : :onresolve
        dirty[] = true
    elseif e.key === Keyboard.space
        paused[] = !paused[]
        dirty[] = true
    else
        return Consume(false)
    end
    return Consume(true)
end

# --- scheduler --------------------------------------------------------------

# The whole non-blocking story. Work is chopped into slices of at most one
# round by construction, with no locks, no channels and no thread-safety
# questions, so the app cannot hang however large the pending set gets.
on(events(fig).tick) do _
    resolved = 0
    if !paused[] && npending(set) > 0
        resolved = run_for!(set, manipulating() ? BUSY_BUDGET : IDLE_BUDGET;
                            chunk = CHUNK)
    end

    # Continuous refinement. Only when everything has resolved: refining on
    # half-finished data wastes the pass, and it would otherwise pile pending
    # samples up faster than the scheduler retires them.
    #
    # Unlike the scheduler, a refinement pass cannot be chopped up — it is one
    # triangulation — and it gets slower as the set grows (~5 ms at 2k
    # samples, ~200 ms at 20k). So rather than holding REFINE_HZ regardless,
    # time the pass and keep the *next* one far enough out that refinement
    # never takes more than REFINE_DUTY of the wall clock. Refinement then
    # slows down by itself as it gets expensive, and the window stays
    # responsive at any sample count. The delay is clamped at REFINE_WAIT
    # because otherwise one freak pass — the first one, which pays for
    # compiling `refine!`, took two seconds — schedules the next one eight
    # seconds out and the loop looks broken.
    #
    # It also yields entirely while the view is being manipulated. A pass
    # cannot be interrupted once started, so at 17k samples an unlucky scroll
    # waits ~200 ms behind one; stepping aside costs nothing and keeps pan and
    # zoom immediate.
    if auto[] && !paused[] && !manipulating() && npending(set) == 0 &&
            length(set) < MAX_SAMPLES && time() >= next_refine[]
        t0 = time()
        refine_view!()
        cost = time() - t0
        next_refine[] = time() + clamp(cost / REFINE_DUTY, 1 / REFINE_HZ, REFINE_WAIT)
    end

    if (resolved > 0 || dirty[]) && time() >= next_commit[]
        commit!(vb[])
        dirty[] = false
        next_commit[] = time() + 1 / COMMIT_HZ
        status[] = string(
            length(set), " samples · ", nresolved(set), " resolved · ",
            npending(set), " pending · ", recon[], " · admit=", vb[].admit,
            auto[] ? " · auto-refine $(REFINE_HZ)Hz" : "",
            length(set) >= MAX_SAMPLES ? " (capped)" : "",
            paused[] ? " · PAUSED" : "")
    end
end

# --- warm-up ----------------------------------------------------------------

"""
Compile every path a click can take, before the user takes one.

Without this the *first* click costs about 15 seconds and the second still
costs two, because Julia is compiling DelaunayTriangulation, Makie's
`poly`/`voronoiplot` pipeline and the GL upload on demand — and it reads as
the app being hung, since nothing is drawn until it finishes. One click's
worth of each distinct path up front takes it to ~300 ms for the first click
and ~20 ms after.

Creating the plot, *updating the generators*, updating only the colours and
running a refinement are separate code paths and each has to be touched;
warming only the first still left the second click at two seconds.
"""
function warmup!()
    corner = (-1.95, -1.85), (0.9, 1.0)      # off in the far exterior
    sprinkle!(set, corner..., 12)
    resolve_all!(set); commit!(vb[]); Makie.colorbuffer(screen)   # create
    sprinkle!(set, corner..., 12; offset = 500)
    commit!(vb[]); Makie.colorbuffer(screen)                      # add generators
    resolve_all!(set); commit!(vb[]); Makie.colorbuffer(screen)   # colours only
    refine!(set, 5; mindist = 1.0e-6)
    resolve_all!(set); commit!(vb[]); Makie.colorbuffer(screen)   # refine + prune
    swap_reconstruction!(); commit!(vb[]); Makie.colorbuffer(screen)   # the other one
    swap_reconstruction!(); commit!(vb[]); Makie.colorbuffer(screen)
    empty!(vb[])
    return nothing
end

screen = display(fig)

if WARMUP
    status[] = "warming up — compiling the tessellation and draw path…"
    sleep(0.1)          # let that frame out before we block the main thread
    @info "scatterbrot: warming up (one-off, ~15 s)"
    @info "scatterbrot: ready in $(round(@elapsed(warmup!()), digits = 1)) s"
    status[] = "0 samples"
end

isinteractive() || wait(screen)
