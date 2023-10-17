export Observables, SentinelObservable

module Observables
    include("abstract_observable.jl")
    include("sentinel_surveillance.jl")
    include("influence_maximization.jl")
    include("immunization.jl")
end


