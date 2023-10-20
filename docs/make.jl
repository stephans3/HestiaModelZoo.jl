using HestiaModelZoo
using Documenter

DocMeta.setdocmeta!(HestiaModelZoo, :DocTestSetup, :(using HestiaModelZoo); recursive=true)

makedocs(;
    modules=[HestiaModelZoo],
    authors="Stephan Scholz",
    repo="https://github.com/stephans3/HestiaModelZoo.jl/blob/{commit}{path}#{line}",
    sitename="HestiaModelZoo.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://stephans3.github.io/HestiaModelZoo.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/stephans3/HestiaModelZoo.jl",
    devbranch="main",
)
