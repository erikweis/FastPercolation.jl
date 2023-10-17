export micro_to_canonical

function binomial_pw(p,N)
    pdf.(Binomial(N,p),1:N)
end

function micro_to_canonical(p::Float64, Qs; M = nothing, null_value = 0)

    M = isa(M,Nothing) ? length(Qs) : M
    b = Binomial(N,p)
    pw = pdf.(b,1:M)
    null_value *= pdf(b,0) # Q(M=0) *pw[0]
    return sum(pw .* Qs) .+ null_value
end

function micro_to_canonical(ps::Vector{Float64}, Qs)
    return map(p->micro_to_canonical(p,Qs),ps)
end