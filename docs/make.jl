using ChromDistABC
using Documenter

DocMeta.setdocmeta!(ChromDistABC, :DocTestSetup, :(using ChromDistABC); recursive=true)

makedocs(;
    modules=[ChromDistABC],
    authors="Hannes Becher",
    repo="https://github.com/hannesbecher/ChromDistABC.jl/blob/{commit}{path}#{line}",
    sitename="ChromDistABC.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
