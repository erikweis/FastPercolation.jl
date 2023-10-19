include("abstract_observable.jl")
include("../utils.jl")
using Graphs

export SentinelObservable, reset_observable, update_observable_data, is_immune, calculate_observables!

QSentinel = Matrix{Int64}

mutable struct SentinelObservable <: Observable
    sentinel_sets::Vector{Vector{Int64}}
    Is::Vector{Vector{Int64}}
    Qs::Vector{QSentinel}
    g::SimpleGraph{Int64}
end

function SentinelObservable(sentinel_sets,n::Int64,m::Int64)
    Is = [zeros(Int64,n) .+ (n+1) for _ in 1:n]
    for i in 1:n; Is[i][i] = 0; end
    Qs = [zeros(Int64,n,length(sentinel_sets)) for _ in 1:m]
    g = SimpleGraph{Int64}(n,0)
    SentinelObservable(sentinel_sets,Is,Qs,g)
end

function reset_observable(M::SentinelObservable)
    n = length(M.Is)
    m = length(M.Qs)
    foreach(1:n, M.Is) do i,I
        I .= (n+1)
        I[i] = 0
    end
    foreach(Q-> Q .= n,M.Qs)
    M.g = SimpleGraph{Int64}(n,0)
end


####### update observable data
function bfs_improve!(root,I,g)
    queue = [root]
    while length(queue)>0
        i = popfirst!(queue)
        for j ∈ neighbors(g,i)
            if (I[j] > I[i] + 1)
                I[j] = I[i] + 1    # node j improves
                push!(queue, j)
            end
        end
    end
end

function update_I!(i, j, I::Vector{Int64},g)
    #check if anything changes
    if I[i] < I[j] - 1      # node i improves its infection time
        bfs_improve!(i,I,g)
    elseif I[j] < I[i] - 1  # node j improves its infection time
        bfs_improve!(j,I,g)
    else
        return nothing      # the infection times are equal or within 1 step, no improvement
    end
end

function update_observable_data(i,j,x::Vector{Int64},M::SentinelObservable)
    add_edge!(M.g,i,j)
    foreach(I -> update_I!(i,j,I,M.g), M.Is)
end

function calculate_observables!(
    m::Int64,
    M::SentinelObservable,
    x::Vector{Int64}
)

    """
    calculate the infected time t over a matrix, where Q[i,j] indicates
    that under seed node i, the set j was infected at time t 
    """
    
    n = length(x)
    Q = M.Qs[m]
    #Q .= n      # set all values to maximum possible, -> ensure reset is run!
    
    for seed in 1:n
        
        I = M.Is[seed]
        max_t = maximum(t for t in I if t < n+1)
        r_seed = find(x,seed)
        for (i, S) in enumerate(M.sentinel_sets)
            for s in S
                if 0 <= I[s] < Q[seed,i] && find(x,s) == r_seed    # check if infection time of s is better, seed must be in same percolation cluster as sentinel
                    Q[seed,i] = I[s]        # node s improves first infection time for S
                end
            end
            if Q[seed,i] > max_t; Q[seed,i] = max_t end    # sentinel set was not infected
        end  
    end
end

is_immune(i,j,M::SentinelObservable) = false