include("abstract_observable.jl")
include("../utils.jl")

export UtilityTransformationObservable, reset_observable, update_observable_data, is_immune, calculate_observables!


mutable struct UtilityTransformationObservable{T} <: Observable
    M::Observable
    Qs::Vector{T}
    Us::Vector{Function}
end

function UtilityTransformationObservable(M::Observable, Us::Vector{Function})
    Qs = [deepcopy(M.Qs) for _ in Us]
    T = typeof(M.Qs)
    UtilityTransformationObservable{T}(M,Qs,Us)
end

reset_observable(M::UtilityTransformationObservable) = nothing

update_observable_data(i,j,x,M::UtilityTransformationObservable) = update_observable_data(i,j,x,M)

function calculate_observables!(
    m::Int64,
    UTM::UtilityTransformationObservable,
    x::Vector{Int64}
)
    calculate_observables!(m,UTM.M,x)
    for (Q, U) in zip(M.Qs, M.Us)
        Q .= U.(UTM.M.Qs[m])
    end
end