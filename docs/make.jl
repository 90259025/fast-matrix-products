include("../src/module.jl")

using Documenter, Main.OurModuleName

makedocs(
    modules = [Main.OurModuleName],
    sitename = "functions.jl",
    doctest = true
)