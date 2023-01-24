function factorial(n::BigInt)
    if n == 0
        return 1
    end
    
    return n * factorial(n - 1)
end

factorial(n::T where {T <: Integer}) = factorial(convert(BigInt, n))


# IN: [a_1,...,a_2n], [b_1,...,b_n]
# OUT: [c_n-1,..., c_2n-1] where c_k = \sum_{i+j=j}a_i b_j
function RMP(A::Vector{<: Union{Integer, Rational{<: Integer}}}, B::Vector{<: Union{Integer, Rational{<: Integer}}})
    n = length(A) >> 1  # n = len(A)/2
    n2 = n >> 1         # n2 = n/2
    
    if n == 1
        return [A[1] * B[1]]
    elseif n2*2 != n    # n is odd
        return 1        #TODO
    end

    α = RMP(A[1:n] .+ A[(n2 + 1):(3*n2)], B[(n2 + 1):end])
    β = RMP(A[(n2 + 1):(3*n2)], B[1:(n2)] .- B[(n2 + 1):end])
    γ = RMP(A[(n2 + 1):(3*n2)] .+ A[(n + 1):2*n], B[1:(n2)])

    return [α + β; γ - β]
end
