export find

function find(x::Vector{Int64}, i::Int64)
    if x[i]<0; return i; end
    return x[i] = find(x,x[i])
end