module BrotViz

using Mandelbrot
using Makie
using ColorSchemes, Colors

include("renderfractal.jl")
include("samples.jl")
include("treecolors.jl")
include("treelayouts.jl")
include("showtree.jl")
include("showspider.jl")
include("showrays.jl")
include("mandelbrotset.jl")
include("juliaset.jl")
include("interiorbinarydecomp.jl")
include("recipes.jl")
include("wedgediagrams.jl")
include("reconstruction.jl")
include("voronoibrot.jl")
include("meshbrot.jl")
include("refine.jl")

export TreeLayout, generation_layout, reingold_tilford_layout,
       embeddedtreeplot, embeddedtreeplot!,
       generationtreeplot, generationtreeplot!,
       treeplot, treeplot!,
       mandelbrotsetplot, mandelbrotsetplot!,
       juliasetplot, juliasetplot!,
       spiderplot!, showspider,
       plotrays, plotrays!,
       dynamicraysplot, dynamicraysplot!,
       kneadingtable,
       expansionwedgeplot, expansionwedgeplot!,
       kneadingwedgeplot, kneadingwedgeplot!,
       expansiondigits, kneadingsymbols,
       SampleStatus, Pending, Escaped, Interior,
       Sample, SampleSet, push_sample!, advance!, run_round!, run_for!,
       resolve_all!, sprinkle!, sprinkle_disc!, npending, nresolved,
       VoronoiBrot, voronoibrotplot, voronoibrotplot!, commit!,
       refine!, banddisagreement, rawdisagreement,
       smoothcount, MeshBrot, meshbrotplot, meshbrotplot!, detach!

end
