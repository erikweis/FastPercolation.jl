include("abstract_observable.jl")
include("../utils.jl")

export UtilityTransformationObservable, reset_observable, update_observable_data, is_immune, calculate_observables!


function convert_type(T::Type{Vector{Vector{Int64}}})
    Vector{Vector{Float64}}
end
function convert_type(T::Type{Vector{Matrix{Int64}}})
    Vector{Matrix{Float64}}
end

mutable struct UtilityTransformationObservable{T} <: Observable
    M::Observable
    Qs::Vector{T}
    Us::Vector{Function}
end

function UtilityTransformationObservable(M::Observable, Us::Vector{Function})
    Qs = [deepcopy(M.Qs) for _ in Us]
    T = convert_type(typeof(M.Qs))

    UtilityTransformationObservable{T}(M,Qs,Us)
end

reset_observable(M::UtilityTransformationObservable) = nothing


update_observable_data(i,j,x,UTM::UtilityTransformationObservable) = update_observable_data(i,j,x,UTM.M)

function calculate_observables!(
    m::Int64,
    UTM::UtilityTransformationObservable,
    x::Vector{Int64}
)
    calculate_observables!(m,UTM.M,x)
    for (Qs, U) in zip(UTM.Qs, UTM.Us)
        Qs[m] .= U.(UTM.M.Qs[m])
    end
end