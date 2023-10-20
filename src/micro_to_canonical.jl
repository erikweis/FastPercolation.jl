import Distributions.Binomial
import Distributions.pdf

export micro_to_canonical, reduce_sentinel_Qs, reduce_sentinel_Q

function binomial_pw(p,N)
    pdf.(Binomial(N,p),1:N)
end

function micro_to_canonical(p::Float64, Qs; M = nothing, null_value = 0)

    M = isa(M,Nothing) ? length(Qs) : M
    b = Binomial(M,p)
    pw = pdf.(b,1:M)
    null_value *= pdf(b,0) # Q(M=0) *pw[0]
    return sum(pw .* Qs) .+ null_value
end

function micro_to_canonical(ps::Vector{Float64}, Qs)
    return map(p->micro_to_canonical(p,Qs),ps)
end


function reduce_sentinel_Qs(Qs)
    map(Qs) do Q
        n, m = size(Q)
        [sum(Q[:,i])/n for i in 1:m]
    end
end

function reduce_sentinel_Q(Q::Matrix{Float64})
    n, m = size(Q)
    [sum(Q[:,i])/n for i in 1:m]
end