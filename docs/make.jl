using StatisticalTests
using Documenter

makedocs(;
    modules=[StatisticalTests],
    authors="Kiwamu Ishikura <ishikura.kiwamu@gmail.com> and contributors",
    repo="https://github.com/i-kiwamu/StatisticalTests.jl/blob/{commit}{path}#L{line}",
    sitename="StatisticalTests.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://i-kiwamu.github.io/StatisticalTests.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/i-kiwamu/StatisticalTests.jl",
)
