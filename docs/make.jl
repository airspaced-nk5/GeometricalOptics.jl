using Documenter
using GeometricalOptics

makedocs(
    sitename = "GeometricalOptics.jl",
    format = Documenter.HTML(),
    modules = [GeometricalOptics],
    pages=[ "index.md",
            "Quickstart.md",
            "Options.md",
            "Examples.md",
            "TroubleConv.md",
            "Reference.md"]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
