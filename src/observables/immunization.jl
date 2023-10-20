include("abstract_observable.jl")
include("../utils.jl")

export ImmunizationObservable, reset_observable, update_observable_data, is_immune, calculate_observables!

QImmunization = Vector{Int64}

mutable struct ImmunizationObservable <: Observable
    immunized_set::Vector{Int64}
    Qs::Vector{QImmunization}
end

function ImmunizationObservable(immunized_set, n,m)
    Qs = [zeros(Int64,n) for _ in 1:m]
    ImmunizationObservable(immunized_set,Qs)
end

reset_observable(M::ImmunizationObservable) = nothing

is_immune(i,j,M::ImmunizationObservable) = (i ∈ M.immunized_set || j ∈ M.immunized_set)

update_observable_data(i,j,x,M::ImmunizationObservable) = nothing

function calculate_observables!(
    m::Int64,
    M::ImmunizationObservable,
    x::Vector{Int64}
)
    """
    calculate the outbreak size for all single-node seeds, for targetted vaccination

    Q is a vector of length one, where Q[i] is the negative outbreak size under seed i
    """

    Q = M.Qs[m]
    for i in 1:length(x)
        Q[i] = x[find(x,i)] # find(i) is the outbreak size of the root of i
    end
end