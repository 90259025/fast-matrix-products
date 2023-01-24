function factorial(n::BigInt)
    if n == 0
        return 1
    end
    
    return n * factorial(n - 1)
end

factorial(n::T) where T <: Integer = factorial(convert(BigInt, n))



# Karastuba style, from https://hal.inria.fr/inria-00071921/document

# IN: [y_1,...,y_n], [z_1,...,z_{2n-1}]
# OUT: [c_{n-1},..., c_{2n-1}] where c_k = \sum_{i + j = k} a_i b_j
function MP(y::Vector{<: Union{Integer, Rational}}, z::Vector{<: Union{Integer, Rational}})
    n = length(y)
    
    if n == 1
        return [y[1] * z[1]]
    end
    
    n0 = n >> 1
    n1 = (n + 1) >> 1
    
    α = MP(y[(n0 + 1):end], z[1:(2n1 - 1)] .+ z[(n1 + 1):(3n1 - 1)])
    
    if iseven(n)
        β = MP(y[(n1 + 1):n] .- y[1:n0], z[(n1 + 1):(3n1 - 1)])
    else
        β = MP([y[n0 + 1]; y[(n1 + 1):n] .- y[1:n0]],z[(n1 + 1):(3n1 - 1)])
    end
    
    γ = MP(y[1:n0], z[(n1 + 1):(n1 + 2n0 - 1)] .+ z[(2n1 + 1):(2n - 1)])
    
    return [α[1:n1] .- β[1:n1]; γ[1:n0] .+ β[1:n0]]
end


# Both delta and Delta from lemma 2 of https://specfun.inria.fr/bostan/publications/BoGaSc07.pdf
# IN: integer d
# OUT: [1/delta(0,d), 1/delta(1,d), ..., 1/delta(d,d)]
function compute_delta(d::T) where T <: Integer
    delta = Array{Rational}(undef, d + 1)
    
    if iseven(d)
        delta[1] = 1 // factorial(d) 
    else
        delta[1] = -1 // factorial(d)
    end
    
    for i = 1:d
        delta[i + 1] = (i - d - 1) // i * delta[i]
    end
    
    return delta
end

# IN: integer a, integer d with a < d
# OUT: [Delta(a,i,d) for i = 0..d]
function compute_Delta(a::T, d::S) where {T <: Integer, S <: Integer}
    Delta = Array{Integer}(undef, d + 1)
    Delta0 = 1
    
    for j = 0:d
        Delta0 *= a - j
    end
    
    Delta[1] = Delta0
    
    for i = 1:d
        Delta[i + 1] = ((a + i)Delta[i]) ÷ (a + i - d - 1)
    end
    
    return Delta
end

# IN: P(0), P(1), ..., P(d); integer a
# OUT: P(a), P(a+1), ..., P(a+d)
# Currrently needs a != d+1, can add that later
function lagrange(p::Vector{<: Union{Integer, Rational}}, a::T) where {T <: Integer}
    d = length(p) - 1

    # obviously cache this at some point
    delta = compute_delta(d)
    Delta = compute_Delta(a,d)

    S = fill(0 // 1, 2d + 1) #padding S... can I save a zero?
    
    p .*= delta

    for i = 0:2d
        S[i + 1] = 1 // (a + i - d)
    end

    prod = MP(p, S) .* Delta
    
    return prod
end
