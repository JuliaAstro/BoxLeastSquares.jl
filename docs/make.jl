using BoxLeastSquares
using Documenter

DocMeta.setdocmeta!(BoxLeastSquares, :DocTestSetup, :(using BoxLeastSquares); recursive=true)

makedocs(;
    modules=[BoxLeastSquares],
    authors="Miles Lucas <mdlucas@hawaii.edu> and contributors",
    repo="https://github.com/mileslucas/BoxLeastSquares.jl/blob/{commit}{path}#{line}",
    sitename="BoxLeastSquares.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mileslucas.github.io/BoxLeastSquares.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mileslucas/BoxLeastSquares.jl",
    devbranch="main",
)
