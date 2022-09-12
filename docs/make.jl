using Documenter
using RoTools

makedocs(
    sitename = "RoTools",
    format = Documenter.HTML(),
    modules = [RoTools]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
