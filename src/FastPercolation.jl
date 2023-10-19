module FastPercolation

#greet() = print("Hello World!")
include("observables/Observables.jl")
#include("observables/sentinel_surveillance.jl")
include("percolation_MC.jl")

export InfluenceMaxObservable, UtilityTransformationObservable, percolation_MC

end # module FastPercolation

