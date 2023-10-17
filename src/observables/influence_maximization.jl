include("abstract_observable.jl")
include("../utils.jl")

export InfluenceMaxObservable, reset_observable, update_observable_data, is_immune, calculate_observables!

QInfluenceMax = Vector{Int64}

mutable struct InfluenceMaxObservable <: Observable
    seed_sets::Vector{Vector{Int64}}
    Qs::Vector{QInfluenceMax}
    roots::Vector{Int64}
end

function InfluenceMaxObservable(seed_sets,n::Int64, m::Int64)
    Qs = [zeros(Int64,length(seed_sets)) for _ in 1:m] 
    max_length = maximum(length(s) for s in seed_sets)
    roots = Vector{Int64}(undef, max_length)
    InfluenceMaxObservable(seed_sets,Qs, roots)
end

reset_observable(M::InfluenceMaxObservable) = nothing

update_observable_data(i,j,x,M::InfluenceMaxObservable) = nothing

is_immune(i,j,M::InfluenceMaxObservable) = false


function calculate_observables!(
    Q::QInfluenceMax,
    M::InfluenceMaxObservable,
    x::Vector{Int64}
)
    """
    calculate the outbreak size rate for each seed set, given seed sets {S}
    Q[i] represents the outbreak under seed set S_i
    """
    
    group = 1
    num_unique = 0

    for (i,S) in enumerate(M.seed_sets)
        num_unique = 0
        for S_j ∈ S
            group = find(x,S_j)
            if group ∉ @view M.roots[1:num_unique]
                num_unique += 1
                M.roots[num_unique] = group
            end
        end
        Q[i] = -sum(@view x[@view M.roots[1:num_unique]])
        #Q[i] = -sum(x[M.roots[1:num_unique]])
    end
    # for (i,S) in enumerate(M.seed_sets)
    #     roots = map(s->find(x,s),S)
    #     unique!(roots)
    #     Q[i] = -sum(x[r] for r in roots)
    # end
end