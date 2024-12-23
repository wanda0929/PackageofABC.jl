using PackageofABC
using Documenter

DocMeta.setdocmeta!(PackageofABC, :DocTestSetup, :(using PackageofABC); recursive=true)

makedocs(;
    modules=[PackageofABC],
    authors="wanda0929",
    sitename="PackageofABC.jl",
    format=Documenter.HTML(;
        canonical="https://wanda0929.github.io/PackageofABC.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/wanda0929/PackageofABC.jl",
    devbranch="main",
)
