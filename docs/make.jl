using BrotViz
using Documenter

makedocs(;
    modules=[BrotViz],
    authors="Jeffrey Wack <jeffwack111@gmail.com> and contributors",
    sitename="BrotViz.jl",
    format=Documenter.HTML(;
        canonical="https://jeffwack.github.io/BrotViz.jl",
        edit_link="main",
        assets=String[],
        size_threshold=nothing,
    ),
    checkdocs=:none,
    pages=[
        "Home" => "index.md",
        "Mandelbrot and BrotViz development notes" => "notes.md",
    ],
)
