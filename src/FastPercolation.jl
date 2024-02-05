module FastPercolation

#greet() = print("Hello World!")
include("observables/Observables.jl")
#include("observables/sentinel_surveillance.jl")
include("percolation_MC.jl")
include("micro_to_canonical.jl")

export InfluenceMaxObservable, UtilityTransformationObservable, percolation_MC, percolation_MC_parallel, micro_to_canonical

end # module FastPercolation
