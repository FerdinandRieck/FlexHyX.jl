using FlexHyX
using Documenter

DocMeta.setdocmeta!(FlexHyX, :DocTestSetup, :(using FlexHyX); recursive=true)

makedocs(;
    modules=[FlexHyX],
    authors="Ferdinand Rieck <ferdinand.rieck@smail.emt.h-brs.de>",
    repo="https://github.com/FerdinandRieck/FlexHyX.jl/blob/{commit}{path}#{line}",
    sitename="FlexHyX.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
