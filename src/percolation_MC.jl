using ..Observables
include("utils.jl")

using Graphs
using ProgressMeter
using Random

export percolation_MC

get_ij(e::Any) = src(e),dst(e)
get_ij(e::Union{Tuple,Vector}) = e[1],e[2]

function percolation_MC!(
    edgelist, 
    x::Vector{Int64},
    M::Observable,
    m_max::Int64
)

    # reset data structures
    reset_observable(M)
    x .= -1

    # shuffle percolation order
    shuffle!(edgelist)
    
    for (m,e) in enumerate(edgelist[1:m_max])

        # get nodes of edge and check for vaccinated nodes
        i,j = get_ij(e)

        #get roots of edges
        r_i = find(x,i)
        r_j = find(x,j)

        # perform union if necessary
        if (r_i !== r_j)
            if x[r_i] < x[r_j]      # size of i is greater than j
                x[r_i] += x[r_j]    # size of i grow by the size of j 
                x[r_j] = r_i        # assign all js to group of i
            else                    # size of i is less than j
                x[r_j] += x[r_i]    # size of j grow by the size of i
                x[r_i] = r_j        # assign all i's to group of j 
            end
        end

        ##### maintain additional data structures necessary for checking qualities
        update_observable_data(i,j,x,M)

        ##### check all qualitiies of interest #######
        calculate_observables!(m,M,x)
    end
end


function percolation_MC(
    g::Graph,
    M::Observable;
    num_samples::Int64 = 10,
    m_max = nothing,
    verbose=false
)

    if m_max === nothing; m_max=ne(g); end
    n = nv(g)
    x = zeros(Int64, n) .- 1 
    edgelist = collect(edges(g))

    # create collector state
    Qs = deepcopy(M.Qs)

    prog = Progress(num_samples)
    for i in 1:num_samples
        percolation_MC!(edgelist, x, M, m_max)
        Qs .+= M.Qs

        if verbose; update!(prog,i);end
    end

    microcanonical = map(Q -> convert.(Float64,Q),Qs)
    return map(Q -> Q./= num_samples, microcanonical)
end


# update_observable_data(i,j,x,M::NeighborhoodObservable) = nothing
# update_observable_data(i,j,x,M::MarginalsObservable) = nothing

# is_immune(i,j,M::MarginalsObservable) = (i ∈ M.immunized_set || j ∈ M.immunized_set)
# is_immune(i,j,M::NeighborhoodObservable) = false
