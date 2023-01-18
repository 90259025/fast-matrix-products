function factorial(n::BigInt)
    if n == 0
        return 1
    end
    
    return n * factorial(n - 1)
end

factorial(n::T where {T <: Integer}) = factorial(convert(BigInt, n))



function RMP(A::Vector{T}, B::Vector{T}) where {T <: Integer}
    n = length(A) ÷ 2
    
    if n == 1
        return [A[1] * B[1]]
    end

    α = RMP(A[1:n] .+ A[(n ÷ 2 + 1):(3n ÷ 2)], B[(n ÷ 2 + 1):end])
    β = RMP(A[(n ÷ 2 + 1):(3n ÷ 2)], B[1:(n ÷ 2)] .- B[(n ÷ 2 + 1):end])
    γ = RMP(A[(n ÷ 2 + 1):(3n ÷ 2)] .+ A[(n + 1):end], B[1:(n ÷ 2)])

    return [α + β; γ - β]
end

A = [1, 2, 3, 4, 8, 9, 26, 7]
B = [4, 6, 38, -1]

println(RMP(A, B))