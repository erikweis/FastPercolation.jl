using Pkg
Pkg.activate(".")
Pkg.add("ParallelUtilities")

#Pkg.instantiate()
Pkg.resolve()


include("test_percolation_MC.jl")