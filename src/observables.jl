# using DrWatson
# @quickactivate "uncertain-influence-max"

# using Graphs
# using Random
# using Dates
# using ProgressMeter
# using Distributions

# abstract type Observable end

# QSentinel = Matrix{Int64}
# QImmunization = Vector{Int64}
# QInfluenceMax = Vector{Int64}
# QNeighborhood = Matrix{Bool}
# QMarginals = Vector{Int64}


# ####### marginals #######
# mutable struct MarginalsObservable <: Observable
#     seed_set::Vector{Int64}
#     immunized_set::Vector{Int64}
#     Qs::Vector{Vector{Int64}}
#     roots::Vector{Int64}
# end
# function MarginalsObservable(seed_set,immunized_set,n::Int64, m::Int64)
#     Qs = [zeros(Int64,n) for _ in 1:m]
#     roots = zeros(Int64,length(seed_set))
#     MarginalsObservable(seed_set,immunized_set,Qs,roots)
# end
# reset_observable(M::MarginalsObservable) = nothing

# ##### three importance problems together
# get_importance_observable(IM::InfluenceMaximization,seed_sets,n,m) = InfluenceMaxObservable(seed_sets,n,m)
# get_importance_observable(OM::OutbreakMinimization,seed_set,n,m) = ImmunizationObservable(seed_set,n,m)
# get_importance_observable(SI::SentinelImportance,seed_sets,n,m) = SentinelObservable(seed_sets,n,m)

# ##### MPL ##########
# mutable struct NeighborhoodObservable <: Observable
#     i::Int64
#     i_idx::Int64
#     neighbors_idx::Vector{Int64}
#     neighborhood_edgelist_idx::Vector{Vector{Int64}}
#     num_nodes::Int64
#     node_map::Vector{Int64}
#     Qs::Vector{Matrix{Bool}}
#     t::Int64
# end
# function NeighborhoodObservable(neighborhood_edgelist, i, num_samples)
    
#     # map edgelist and node to an index
#     node_map = sort(unique(vcat(neighborhood_edgelist...)))
#     neighbor_nodes = setdiff(node_map,i)

#     num_nodes = length(node_map)
#     num_neighbor_nodes = length(neighbor_nodes)
#     f(x,node_map) = findfirst(y->y==x,node_map)
#     i_idx = f(i,node_map)
#     #println("num_nodes: $num_nodes, i_idx: $i_idx")
#     neighbors_idx = (i_idx === nothing) ? (1:num_nodes) : setdiff(1:num_nodes,i_idx)
#     neighborhood_edgelist_idx = [map(x->f(x,node_map),ij) for ij in neighborhood_edgelist]
#     # define observable
#     m = length(neighborhood_edgelist)
#     Qs = [zeros(Bool,(num_samples,num_neighbor_nodes)) for _ in 1:m] #Qs[m][t,j] is a boolean, storing 1 if node j is reachable from node i at the t-th sample with m edges placed

#     NeighborhoodObservable(i,i_idx,neighbors_idx, neighborhood_edgelist_idx,num_nodes,node_map,Qs,1)
# end
# reset_observable(M::NeighborhoodObservable) = nothing
# # function reset_observable(M::NeighborhoodObservable)
# #     remove_edges!(M.g, edges(M.g))
# # end

# mutable struct CompoundObservable <: Observable
#     Qs::Vector{Observable}
# end

# function find(x::Vector{Int64}, i::Int64)
#     if x[i]<0; return i; end
#     return x[i] = find(x,x[i])
# end

# function bfs_improve!(root,I,g)
#     queue = [root]
#     while length(queue)>0
#         i = popfirst!(queue)
#         for j ∈ neighbors(g,i)
#             if (I[j] > I[i] + 1)
#                 I[j] = I[i] + 1    # node j improves
#                 push!(queue, j)
#             end
#         end
#     end
# end


# function update_I!(i, j, I::Vector{Int64},g)
#     #check if anything changes
#     if I[i] < I[j] - 1      # node i improves its infection time
#         bfs_improve!(i,I,g)
#     elseif I[j] < I[i] - 1  # node j improves its infection time
#         bfs_improve!(j,I,g)
#     else
#         return nothing      # the infection times are equal or within 1 step, no improvement
#     end
# end



# update_observable_data(i,j,x,M::NeighborhoodObservable) = nothing
# update_observable_data(i,j,x,M::ImmunizationObservable) = nothing
# update_observable_data(i,j,x,M::InfluenceMaxObservable) = nothing
# update_observable_data(i,j,x,M::MarginalsObservable) = nothing

# is_immune(i,j,M::ImmunizationObservable) = (i ∈ M.immunized_set || j ∈ M.immunized_set)
# is_immune(i,j,M::MarginalsObservable) = (i ∈ M.immunized_set || j ∈ M.immunized_set)
# is_immune(i,j,M::InfluenceMaxObservable) = false
# is_immune(i,j,M::SentinelObservable) = false
# is_immune(i,j,M::NeighborhoodObservable) = false



# function calculate_observables!(
#     Q::QImmunization,
#     M::ImmunizationObservable,
#     x::Vector{Int64}
# )
#     """
#     calculate the outbreak size for all single-node seeds, for targetted vaccination

#     Q is a vector of length one, where Q[i] is the negative outbreak size under seed i
#     """
#     for i in 1:length(x)
#         Q[i] = x[find(x,i)] # find(i) is the outbreak size of the root of i
#     end
# end

# function calculate_observables!(
#     Q::QInfluenceMax,
#     M::InfluenceMaxObservable,
#     x::Vector{Int64}
# )
#     """
#     calculate the outbreak size rate for each seed set, given seed sets {S}
#     Q[i] represents the outbreak under seed set S_i
#     """
    
#     group = 1
#     num_unique = 0

#     for (i,S) in enumerate(M.seed_sets)
#         num_unique = 0
#         for S_j ∈ S
#             group = find(x,S_j)
#             if group ∉ @view M.roots[1:num_unique]
#                 num_unique += 1
#                 M.roots[num_unique] = group
#             end
#         end
#         Q[i] = -sum(@view x[@view M.roots[1:num_unique]])
#         #Q[i] = -sum(x[M.roots[1:num_unique]])
#     end
#     # for (i,S) in enumerate(M.seed_sets)
#     #     roots = map(s->find(x,s),S)
#     #     unique!(roots)
#     #     Q[i] = -sum(x[r] for r in roots)
#     # end
# end

# function calculate_observables!(
#     Q::QMarginals,
#     M::MarginalsObservable,
#     x::Vector{Int64}
# )
#     """
#     calculate whether each node is reachable from the seed set. Since no edge has been added that connects a
#     vaccinated node, this should automatically ensure vaccinated nodes are uninfected.
#     """

#     for (i,s) in enumerate(M.seed_set)
#         M.roots[i] = find(x,s)
#     end
#     #println(M.roots)
#     #println(Q)
#     for i in eachindex(Q)
#         Q[i] = find(x,i) ∈ M.roots
#     end
# end

# function calculate_observables!(
#     Q::QNeighborhood,
#     M::NeighborhoodObservable,
#     x::Vector{Int64}
# )

#     t = M.t
#     foreach(enumerate(M.neighbors_idx)) do (idx,j)
#         Q[t,idx] = find(x,M.i_idx) == find(x,j)
#     end
# end

# get_ij(e::Any) = src(e),dst(e)
# get_ij(e::Union{Tuple,Vector}) = e[1],e[2]

# function percolation_MC!(
#     edgelist, 
#     x::Vector{Int64},
#     M::Observable,
#     m_max::Int64
# )

#     # reset data structures
#     reset_observable(M)
#     x .= -1

#     # shuffle percolation order
#     shuffle!(edgelist)
    
#     for (idx,e) in enumerate(edgelist[1:m_max])

#         # get nodes of edge and check for vaccinated nodes
#         i,j = get_ij(e)
#         if !is_immune(i,j,M)            # edge is active only when neither neither node is immune

#             #get roots of edges
#             r_i = find(x,i)
#             r_j = find(x,j)

#             # perform union if necessary
#             if (r_i !== r_j)
#                 if x[r_i] < x[r_j]      # size of i is greater than j
#                     x[r_i] += x[r_j]    # size of i grow by the size of j 
#                     x[r_j] = r_i        # assign all js to group of i
#                 else                    # size of i is less than j
#                     x[r_j] += x[r_i]    # size of j grow by the size of i
#                     x[r_i] = r_j        # assign all i's to group of j 
#                 end
#             end
#         end

#         ##### maintain additional data structures necessary for checking qualities
#         update_observable_data(i,j,x,M)

#         ##### check all qualitiies of interest #######
#         calculate_observables!(M.Qs[idx],M,x)
#     end
# end


# function percolation_MC(
#     g::Graph,
#     M::Observable;
#     num_samples::Int64 = 10,
#     m_max = nothing,
#     verbose=false
# )
#     if m_max === nothing; m_max=ne(g); end
#     n = nv(g)
#     x = zeros(Int64, n) .- 1 
#     edgelist = collect(edges(g))

#     # create collector state
#     Qs = deepcopy(M.Qs)

#     prog = Progress(num_samples)
#     for i in 1:num_samples
#         percolation_MC!(edgelist, x, M, m_max)
#         Qs .+= M.Qs

#         if verbose; update!(prog,i);end
#     end

#     microcanonical = map(Q -> convert.(Float64,Q),Qs)
#     return map(Q -> Q./= num_samples, microcanonical)
# end



# function binomial_pw(p,N)
#     pdf.(Binomial(N,p),1:N)
#     # pw = zeros(Float64,N)
    
#     # # set un-normalized peak value
#     # n_max = Int(round(p*N))
#     # n_max = maximum([n_max,1])
#     # pw[n_max] = 1
#     # for n in n_max+1:N
#     #     pw[n] = pw[n-1]*((N-n+1)/n)*(p/(1-p))
#     # end
#     # for n in reverse(1:n_max-1)
#     #     pw[n] = pw[n+1]*((n+1)/(N-n))*(1-p)/p
#     # end
#     # return pw ./ (sum(pw) + (1-p)^N )
# end

# function micro_to_canonical(p::Float64, Qs; null_value = 0)

#     N = length(Qs)
#     b = Binomial(N,p)
#     pw = pdf.(b,1:N)
#     null_value *= pdf(b,0) # Q(M=0) *pw[0]
#     return sum(pw .* Qs) .+ null_value
# end

# function micro_to_canonical(ps::Vector{Float64}, Qs)
#     return map(p->micro_to_canonical(p,Qs),ps)
# end


# function reduce_sentinel_Qs(Qs)
#     map(Qs) do Q
#         n, m = size(Q)
#         [sum(Q[:,i])/n for i in 1:m]
#     end
# end

# function reduce_sentinel_Q(Q::Matrix{Float64})
#     n, m = size(Q)
#     [sum(Q[:,i])/n for i in 1:m]
# end
