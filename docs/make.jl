using UncertainHistogramming
using Documenter

DocMeta.setdocmeta!(UncertainHistogramming, :DocTestSetup, :(using UncertainHistogramming); recursive=true)

makedocs(;
    modules=[UncertainHistogramming],
    authors="W. Joe Meese <meese022@umn.edu> and contributors",
    repo="https://github.com/meese-wj/UncertainHistogramming.jl/blob/{commit}{path}#{line}",
    sitename="UncertainHistogramming.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://meese-wj.github.io/UncertainHistogramming.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API Reference" => "api.md"
    ],
)

deploydocs(;
    repo="github.com/meese-wj/UncertainHistogramming.jl",
    devbranch="main",
)
