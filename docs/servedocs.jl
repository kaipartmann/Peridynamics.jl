#-------------------------------------------------------------------------------------------
# RUN THIS SCRIPT ONLY FOR LOCAL DEVELOPMENT OF THE DOCUMENTATION!
#-------------------------------------------------------------------------------------------

# Root of the repository
const REPO_ROOT = dirname(@__DIR__)

# Make sure the docs environment is active and instantiated
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# Communicate with make.jl that docs are build in live mode
push!(ARGS, "LIVE_MODE")

using LiveServer
servedocs(;
    # Documentation root where make.jl and src/ are located
    foldername = joinpath(REPO_ROOT, "docs"),
    # Watch the src folder so docstrings can be Revise'd
    include_dirs = [joinpath(REPO_ROOT, "src")],
    # Skip the folder where Literate.jl output is written. This is needed
    # to avoid infinite loops where running make.jl updates watched files,
    # which then triggers a new run of make.jl etc...
    skip_dirs = [joinpath(REPO_ROOT, "docs", "src", "generated")],
)
