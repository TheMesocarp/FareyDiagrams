using FareyDiagrams 
using Documenter


DocMeta.setdocmeta!(FareyDiagrams, :DocTestSetup, :(using FareyDiagrams); recursive=true)

makedocs(
    modules = [FareyDiagrams], 
    sitename = "FareyDiagrams.jl", 
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://Adalovescoffee.github.io/FareyDiagrams.jl/stable" 
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Introduction" => "man/introduction.md",
            "Functions" => "man/functions.md", 
        ],
        
        "API Reference" => "api.md", 
    ],
   
    warnonly = [:missing_docs, :cross_references],
)


deploydocs(
    repo = "github.com/Adalovescoffee/FareyDiagrams.jl.git", 
    devbranch = "main", 
    push_preview = true, 
)