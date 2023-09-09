try
    @eval using Documenter, MRINavigator
catch e
    Pkg.develop(PackageSpec(path = joinpath(@__DIR__, "..")))
    Pkg.instantiate()
    @eval using Documenter, MRINavigator
end

 
makedocs(modules=[MRINavigator],
        sitename = "MRINavigator.jl",
        authors = "Laura Beghini",
        pages = [
        "Home" => "index.md",
        "Getting Started" => "GettingStarted.md",
        "Pipelines" => "Pipelines.md",
        "API" => "API.md"
        ],
        )
 
deploydocs(;
    repo = "github.com/NordicMRspine/MRINavigator",
    push_preview = true,
    deploy_config = Documenter.GitHubActions(),
    )
