using BoxLeastSquares
using Documenter

DocMeta.setdocmeta!(BoxLeastSquares, :DocTestSetup, :(using BoxLeastSquares); recursive=true)
include("pages.jl")
makedocs(;
    modules=[BoxLeastSquares],
    authors="Miles Lucas <mdlucas@hawaii.edu> and contributors",
    repo="https://github.com/JuliaAstro/BoxLeastSquares.jl/blob/{commit}{path}#{line}",
    sitename="BoxLeastSquares.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaAstro.github.io/BoxLeastSquares.jl",
        assets=String[],
    ),
    pages=pages
)

deploydocs(;
    repo="github.com/JuliaAstro/BoxLeastSquares.jl",
    push_preview=true,
    devbranch="main",
)
