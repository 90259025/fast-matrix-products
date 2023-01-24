function factorial(n::BigInt)
    if n == 0
        return 1
    end
    
    return n * factorial(n - 1)
end

factorial(n::T where {T <: Integer}) = factorial(convert(BigInt, n))



# Karastuba style, from https://hal.inria.fr/inria-00071921/document

# IN: [y_1,...,y_n], [z_1,...,z_2n-1]
# OUT: [c_n-1,..., c_2n-1] where c_k = \sum_{i+j=j}a_i b_j
function MP(y::Vector{<: Union{Integer, Rational{<: Integer}}}, z::Vector{<: Union{Integer, Rational{<: Integer}}})
    n = length(y)
    if n == 1
        return [y[1] * z[1]]
    end
    n0 = n >> 1
    n1 = (n+1) >> 1
    α = MP(y[n0+1:end], z[1:2*n1-1] .+ z[n1+1:3*n1-1])
    if iseven(n)
        β = MP(y[n1+1:n] .- y[1:n0], z[n1+1:3*n1-1])
    else
        β = MP([y[n0+1]; y[n1+1:n] .- y[1:n0]],z[n1+1:3*n1-1])
    end
    γ = MP(y[1:n0], z[n1+1:n1+2*n0 - 1] .+ z[2*n1+1:2*n-1])
    return [α[1:n1] .- β[1:n1]; γ[1:n0] .+ β[1:n0]]
end