using BoxLeastSquares
using Documenter

DocMeta.setdocmeta!(BoxLeastSquares, :DocTestSetup, :(using BoxLeastSquares); recursive=true)
include("pages.jl")
makedocs(;
    modules = [BoxLeastSquares],
    authors = "Miles Lucas <mdlucas@hawaii.edu> and contributors",
    repo = "https://github.com/JuliaAstro/BoxLeastSquares.jl/blob/{commit}{path}#{line}",
    sitename = "BoxLeastSquares.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        assets = String[],
        canonical = "https://JuliaAstro.org/BoxLeastSquares/stable/",
    ),
    pages = pages
)

deploydocs(;
    repo = "github.com/JuliaAstro/BoxLeastSquares.jl",
    push_preview = true,
    devbranch = "main",
    versions = ["stable" => "v^", "v#.#"], # Restrict to minor releases
)
