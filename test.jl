<<<<<<< HEAD
include("functions.jl")

# test RMP
=======
function RMP(A::Vector{<: Union{Integer, Rational{<: Integer}}}, B::Vector{<: Union{Integer, Rational{<: Integer}}})
    n = length(A) ÷ 2
    
    if n == 1
        return [A[1] * B[1]]
    end

    α = RMP(A[1:n] .+ A[(n ÷ 2 + 1):(3n ÷ 2)], B[(n ÷ 2 + 1):end])
    β = RMP(A[(n ÷ 2 + 1):(3n ÷ 2)], B[1:(n ÷ 2)] .- B[(n ÷ 2 + 1):end])
    γ = RMP(A[(n ÷ 2 + 1):(3n ÷ 2)] .+ A[(n + 1):end], B[1:(n ÷ 2)])

    return [α + β; γ - β]
end
>>>>>>> 21d35e1b7418885c6cc7c9ea818c33dbd3ae4738

# 7*t^7 + 26*t^6 + 9*t^5 + 8*t^4 + 4*t^3 + 3*t^2 + 2*t + 1
A = [1, 2, 3, 4, 8, 9, 26, 7]
# -t^3 + 38*t^2 + 6*t + 4
B = [4, 6, 38, -1]
# C = -7*t^10 + 240*t^9 + 1021*t^8 + 518*t^7 + 458*t^6 + 233*t^5 + 168*t^4 + 109*t^3 + 62*t^2 + 14*t + 4
println([109, 168, 233, 458] == RMP(A, B))

# 880*t^4 + 4*t^3 + 3*t^2 + 2*t + 1
# A = [1,2,3,4,880,0,0,0,0,0]
# 9*t^4 + 9*t^3 + 8*t^2 + 3*t + 1
# B = [1,3,8,9,9]
# 7920*t^8 + 7956*t^7 + 7103*t^6 + 2717*t^5 + 943*t^4 + 38*t^3 + 17*t^2 + 5*t + 1
# println(RMP(A,B)==[943, 2717, 7103, 7956, 7920])
println(RMP(A,B))