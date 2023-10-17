
using DrWatson
@quickactivate "uncertain-influence-max"

using Graphs
using Random
using Dates
using ProgressMeter
using Distributions
using Statistics

include(srcdir("contagion_models/percolation_monte_carlo.jl"))


abstract type ConfidenceAnalysis end

mutable struct BestSeparationCA{DiffType} <: ConfidenceAnalysis
    canonical::Vector{Float64}
    pw::Vector{Float64}
    infection_rate::Float64
    diff::DiffType
    var::DiffType
    p_values::Vector{Float64}
    all_p_values::Vector{Vector{Float64}}
    save_data::Bool
end

function BestSeparationCA(m::Int64, infection_rate::Float64, M::InfluenceMaxObservable; save_data = false)
    canonical = Vector{Float64}(undef,length(M.seed_sets))
    pw = binomial_pw(infection_rate,m)
    p_values = deepcopy(canonical)
    all_p_values = []
    diff = Vector{Float64}(undef,length(M.seed_sets))
    var = deepcopy(diff)
    save_data = save_data
    BestSeparationCA{Vector{Float64}}(canonical, pw, infection_rate,diff,var,p_values,all_p_values,save_data)
end

mutable struct TrackVariance{QType} <: ConfidenceAnalysis
    canonical::Vector{Float64}
    pw::Vector{Float64}
    infection_rate::Float64
    expected_mean::QType
    expected_mean2::Vector{Float64}
    expected_var::Vector{Float64}
    total_vars::Vector{Float64}
    means::Vector{Float64}
    all_means::Vector{Vector{Float64}}
    all_SEs::Vector{Vector{Float64}}
    save_data::Bool
end
function TrackVariance(m::Int64, infection_rate::Float64, M::InfluenceMaxObservable; save_data = false)
    canonical = Vector{Float64}(undef,length(M.seed_sets))
    pw = binomial_pw(infection_rate,m)
    # diff = Vector{Float64}(undef,length(M.seed_sets))
    # var = deepcopy(diff)
    # expected_mean = zeros(Float64,length(M.seed_sets))
    # expected_mean2 = zeros(Float64,length(M.seed_sets))
    # expected_var = zeros(Float64,length(M.seed_sets))
    total_vars = zeros(Float64,length(M.seed_sets))
    means = zeros(Float64,length(M.seed_sets))
    SEs = zeros(Float64,length(M.seed_sets))
    all_means = []
    all_SEs = []
    save_data = save_data
    TrackVariance{Vector{Float64}}(canonical, pw, infection_rate,expected_mean,expected_mean2,expected_var,all_means,all_SEs,save_data)
end
function TrackVariance(m::Int64, infection_rate::Float64, M::SentinelObservable; save_data = false)
    canonical = Vector{Float64}(undef,length(M.seed_sets))
    pw = binomial_pw(infection_rate,m)
    diff = Vector{Float64}(undef,length(M.seed_sets))
    var = deepcopy(diff)
    mean_means = zeros(Float64,length(M.seed_sets))
    var_means = zeros(Float64,length(M.seed_sets))
    mean_vars = zeros(Float64,length(M.seed_sets))
    means = zeros(Float64,length(M.seed_sets))
    SEs = zeros(Float64,length(M.seed_sets))
    all_means = []
    all_SEs = []
    save_data = save_data
    TrackVariance{Vector{Float64}}(canonical, pw, infection_rate,mean_means,var_means, mean_vars,SEs,all_means,all_SEs,save_data)
end


function calculate_confidence!(
    CA::BestSeparationCA, Qs_mean, Qs_ME2, tdist, num_samples; 
    verbose = true
)

    # define variables
    diff = CA.diff
    var = CA.var
    pw = CA.pw
    canonical = CA.canonical
    p_values = CA.p_values

    # calculate ranking
    micro_to_canonical!(pw,canonical,Qs_mean)
    best_idx = argmax(canonical)

    # calculate p_values for each step of the algorithm separately
    fill!(p_values,0.0)
    for (m_i,pw_i) in zip(1:m,pw)

        # calculate t-statistic
        calculate_diff!(diff,Qs_mean[m_i],best_idx)
        calculate_variance!(var, Qs_ME2[m_i],best_idx,num_samples)
        diff ./= var
        t_statistic = diff

        # calculate p_value from t-statistic
        ps = cdf.(Ref(tdist), t_statistic)
        ps[best_idx] = 0.0
        ps .*= pw_i
        ps[isnan.(ps)] .= 0 # set small values to zero
        p_values .+= ps 
    end

    # print information if necessary
    if verbose
        second_best_quality = sort(canonical, rev=true)[2]
        best_quality = sort(canonical, rev=true)[1]
        dq = best_quality - second_best_quality
        println("\n $num_samples max p_value: $(maximum(p_values)), best_set:$(M.seed_sets[best_idx]), best_q:$best_quality, diff: $dq")
    end

    # keep track over time if desired
    if CA.save_data
        push!(CA.all_p_values,deepcopy(p_values))
    end

end

function calculate_confidence!(
    CA::TrackVariance, Qs_mean, Qs_ME2, tdist, num_samples; 
    verbose = true
)

    # define variables
    pw = CA.pw
    all_SEs = CA.all_SEs
    all_means = CA.all_means
    n = num_samples
    # expected_mean = CA.expected_mean
    # expected_mean2 = CA.expected_mean2
    # expected_var = CA.expected_var
    total_vars = 
    means = CA.means
    SEs = CA.SEs

    # calculate p_values for each step of the algorithm separately
    expected_mean = zeros(Float64,size(Qs_mean[1]))
    expected_mean2 = deepcopy(expected_mean)
    expected_var = deepcopy(expected_mean)
    # fill!(expected_mean,0.0)
    # fill!(expected_mean2,0.0)
    # fill!(expected_var,0.0)

    foreach(1:m) do i

        mean = Qs_mean[i]
        vars = Qs_ME2[i] ./ (n-1)

        expected_mean .+= pw[i] .* mean
        expected_mean2 .+= pw[i] .* (mean .^ 2)
        expected_var .+= pw[i] .* vars

    end
    total_vars .= expected_var .+ expected_mean2 .- (expected_mean .^ 2)
    SEs .= sqrt.( total_vars  ./ n)

    # print information if necessary
    if verbose
        println("\nSEs: $SEs")
        println("\ntotal variances: $total_vars")
    end

    # keep track over time if desired
    if CA.save_data
        push!(all_means,deepcopy(expected_mean))
        push!(all_SEs,deepcopy(SEs))
    end

end



function percolation_MC_adaptive(
    g::Graph,
    M::Observable;
    confidence = 0.05,
    check_interval = 10^2,
    min_samples = 10^2,
    max_samples = 10^3,
    CAs = nothing
)
    
    n = nv(g)
    m = ne(g)
    x = zeros(Int64, n) .- 1 
    edgelist = collect(edges(g))

    # create collector state
    Qs_mean = map(Q -> convert.(Float64,Q), deepcopy(M.Qs) )
    Qs_ME2 = deepcopy(Qs_mean)
    delta = deepcopy(Qs_mean)
    delta2 = deepcopy(Qs_mean)


    @showprogress for i in 1:max_samples
        
        # do percolation
        percolation_MC!(edgelist, x, M)
        
        # update mean and variance
        delta .= M.Qs .- Qs_mean
        Qs_mean .+= (delta ./ i)
        delta2 .= M.Qs .- Qs_mean
        for (Q_ME2, d1, d2) in zip(Qs_ME2, delta, delta2)
            Q_ME2 .+= d1 .* d2
        end
        
        # check p_values every so often
        if i % check_interval == 0 && (i >= min_samples)
            
            tdist = TDist(2*i-2)
            
            for CA in CAs
                calculate_confidence!(CA, Qs_mean, Qs_ME2, tdist, i)
            end
        end
    end
    Qs_mean
end


function get_best_idx(canonical::Vector{Float64})
    argmax(canonical)
end


function get_best_idx(canonical::Matrix{Float64})
    argmax(sum(canonical,dims=1))[2] # the argmax function returns coordinates of 1xK matrix (1,idx)
end


function calculate_diff!(
    diff::Vector{Float64},
    Q::Vector{Float64},
    best_idx::Int64
)
    diff .= Q .- Q[best_idx]
end

function calculate_SE!(var::Vector{Float64},Q_ME2::Vector{Float64}, n::Int64)
    var .= Q_ME2 ./ (n*(n-1))
    sqrt.(var)
end

function calculate_variance!(
    var::Vector{Float64},
    Q_ME2::Vector{Float64},
    best_idx::Int64,
    n::Int64 # num samples
)
    var .= (Q_ME2 .+ Q_ME2[best_idx])
    var .= sqrt.(var)
    var ./= n
end


function calculate_diff!(
    diff::Matrix{Float64},
    Q::Matrix{Float64},
    best_idx::Int64
)
    diff .= Q .- Q[:,best_idx]
end


function micro_to_canonical!(pw::Vector{Float64},canonical::Vector{Float64},Qs)
    canonical .= sum(pw .* Qs)
end