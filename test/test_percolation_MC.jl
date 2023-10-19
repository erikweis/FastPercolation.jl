
using FastPercolation
using Test
using Graphs
using Random
using StatsBase
using TimerOutputs

function test_percolation_MC()

    n = 20
    g = newman_watts_strogatz(n, 4, 0.3)
    Random.seed!(1)
    num_sets = 50
    rsets = [sample(1:n,3,replace=false) for _ in 1:num_sets]
    M = SentinelObservable(rsets,n,ne(g))
    Us::Vector{Function} = [x->x]
    UtilityTransformationObservable(M,Us)
    #edgelist = collect(edges(g))
    x = zeros(Int64, n) .- 1 
    out = percolation_MC(g,M; num_samples = 10^1)
    print(out[1])
    return true
end

@testset begin
    @test test_percolation_MC()
end

# function test_create_Ms()

#     sets = [[1,2,3],[4,5,6]]
#     M_sent = SentinelObservable(sets,10,20)
#     M_imm = ImmunizationObservable(sets[1],10,20)
#     M_inf = InfluenceMaxObservable(sets,10,20)
#     M_comp = CompoundObservable([M_imm,M_inf])

#     return true

# end

# function test_find()
#     x = [-1,1,2,3,4,4]
#     println(find(x,5))
#     return true
# end


# function basic_graph()
#     edges = [(0,1),(1,2),(2,3),(3,4),(3,5),(4,5)]
#     g = SimpleGraph{Int64}(5,0)
#     for (i,j) in edges
#         add_edge!(g,i,j)
#     end
#     g
# end


# function test_bfs_improve()

#     g = basic_graph()
#     I = zeros(Int,5) .+ 5
#     I[1] = 0
#     bfs_improve!(1,I,g)
#     check1 = I == [0,1,2,3,3]

#     add_edge!(g,5,2)
#     bfs_improve!(2,I,g)
#     check2 = I == [0,1,2,3,2]
#     return check1 && check2

# end

# function test_update_I()

#     g = basic_graph()
#     i,j = 5,1
    
#     I = zeros(Int,5) .+ 5
#     I[1] = 0
#     bfs_improve!(1,I,g)
#     println(I)

#     add_edge!(g,i,j)
#     update_I!(i,j,I,g)
#     println(I)
#     I == [0,1,2,2,1]
# end

# function test_update_observable_data()

#     x = [-1,-1,-1,-1,-1,-1]
#     I[1] = 0
#     bfs_improve!(1,I,g)

#     sets = [[1,2,3],[4,5],[2,3],[1],[2],[3],[4],[5]]
#     M_sent = SentinelObservable(sets,10)
    

#     update_observable_data(1,2)

# end


# function test_percolation_MC!()

#     n = 10
#     g = newman_watts_strogatz(n, 4, 0.3)
#     M = InfluenceMaxObservable([[i] for i in 1:n],n,ne(g))
#     M2 = SentinelObservable([[i] for i in 1:5],n,ne(g))
#     M3 = ImmunizationObservable([1,4,5],n,ne(g))
       
#     edgelist = collect(edges(g))
#     x = zeros(Int64, n) .- 1 
#     #percolation_MC!(edgelist,x,M)
#     percolation_MC!(edgelist,x,M2)
#     #percolation_MC!(edgelist,x,M3)

#     return true
# end




# function test_percolation_MC_sentinel()

#     n = 200
#     g = newman_watts_strogatz(n, 4, 0.3)
#     Random.seed!(1)
#     num_sets = 50
#     rsets = [sample(1:n,20,replace=false) for _ in 1:num_sets]
#     M = SentinelObservable(rsets,n,ne(g))
    
#     edgelist = collect(edges(g))
#     x = zeros(Int64, n) .- 1 
#     percolation_MC(g,M; num_samples = 10^3)
    
#     return true
# end

# function test_micro_to_canonincal()

#     n = 10
#     g = newman_watts_strogatz(n, 4, 0.3)
#     M = SentinelObservable([[i] for i in 1:n],n,ne(g))
    
#     edgelist = collect(edges(g))
#     x = zeros(Int64, n) .- 1 
#     microcanonical = percolation_MC(g,M; num_samples = 10^3)
    
#     ps = collect(0.0:0.05:1)
#     out = micro_to_canonical(ps, microcanonical)
#     println(out)

#     quality_s2 = [ x[1] for x in out ]
#     println(quality_s2)
# end

# # @test test_create_Ms()
# # @test test_find()
# # @test test_bfs_improve()
# # @test test_update_I()
# #@test test_percolation_MC!()

# global to = TimerOutput()
# @timeit to "inf max" test_percolation_MC()
# show(to)

# global to = TimerOutput()
# @timeit to "sentinel importance" test_percolation_MC_sentinel()
# show(to)

