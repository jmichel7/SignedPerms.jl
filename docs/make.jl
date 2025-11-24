using Documenter, SignedPerms, PermGroups

DocMeta.setdocmeta!(SignedPerms, :DocTestSetup, :(using SignedPerms); recursive=true)

makedocs(;
    modules=[SignedPerms],
    authors="Jean Michel <jean.michel@imj-prg.fr> and contributors",
    sitename="SignedPerms.jl",
    format=Documenter.HTML(;
        canonical="https://jmichel7.github.io/SignedPerms.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    warnonly=:missing_docs,
)

deploydocs(;
    repo="github.com/jmichel7/SignedPerms.jl",
    devbranch="main",
)
