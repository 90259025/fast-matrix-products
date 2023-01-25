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
    
    n₀ = n >> 1
    n₁ = (n + 1) >> 1
    
    α = MP(y[(n₀ + 1):end], z[1:(2n₁ - 1)] .+ z[(n₁ + 1):(3n₁ - 1)])
    
    if iseven(n)
        β = MP(y[(n₁ + 1):n] .- y[1:n₀], z[(n₁ + 1):(3n₁ - 1)])
    else
        β = MP([y[n₀ + 1]; y[(n₁ + 1):n] .- y[1:n₀]],z[(n₁ + 1):(3n₁ - 1)])
    end
    
    γ = MP(y[1:n₀], z[(n₁ + 1):(n₁ + 2n₀ - 1)] .+ z[(2n₁ + 1):(2n - 1)])
    
    return [α[1:n₁] .- β[1:n₁]; γ[1:n₀] .+ β[1:n₀]]
end



# Both delta and Delta from lemma 2 of https://specfun.inria.fr/bostan/publications/BoGaSc07.pdf
# IN: integer d
# OUT: [1/delta(0,d), 1/delta(1,d), ..., 1/delta(d,d)]
function compute_δ(d::T) where T <: Integer
    δ = Array{Rational}(undef, d + 1)
    
    if iseven(d)
        δ[1] = 1 // factorial(d) 
    else
        δ[1] = -1 // factorial(d)
    end
    
    for i = 1:d
        δ[i + 1] = (i - d - 1) // i * δ[i]
    end
    
    return δ
end

compute_delta = compute_δ



# IN: integer a, integer d with a < d
# OUT: [Delta(a,i,d) for i = 0..d]
function compute_Δ(a::T, d::S) where {T <: Integer, S <: Integer}
    Δ = Array{Integer}(undef, d + 1)
    Δ₀ = 1
    
    for j = 0:d
        Δ₀ *= a - j
    end
    
    Δ[1] = Δ₀
    
    for i = 1:d
        Δ[i + 1] = ((a + i)Δ[i]) ÷ (a + i - d - 1)
    end
    
    return Δ
end

compute_Delta = compute_Δ


# IN: P(0), P(1), ..., P(d); integer a
# OUT: P(a), P(a+1), ..., P(a+d)
# Only supports a > d, but that suffices
function lagrange(p::Vector{<: Union{Integer, Rational}}, a::T) where {T <: Integer}
    d = length(p) - 1

    if a == d+1
        # use theorem 3.1 from https://dl.acm.org/doi/pdf/10.1145/120694.120697
        G = fill(0//1,2*d+4)
        for i = 1:d+2
            G[i] = (2*(i&1)-1) * binomial(d+1,i-1)
            # G[i] = (-1)^(i-1) * binomial(d+1,i-1)
        end
        H = fill(0//1,2*d+4)
        for i = 1:d+3 #I think 1:d+1 works too but I'm a little nervous about changing it
            H[d+1+i] = -(2*(i&1)-1) * binomial(-d-1,i-1)
        end
        F1 = [0;MP([p;0],G)[1:end-1]]
        return MP(F1,H)[2:end]
    end

    # obviously cache this at some point
    δ = compute_δ(d)
    Δ = compute_Δ(a, d)

    p_tilde = p .* δ

    S = fill(0 // 1, 2d + 1)
    
    for i = 0:2d
        S[i + 1] = 1 // (a + i - d)
    end

    #Q_k is the coefficient of x^{k + d} in ̃pS, and P(a + k) = Q_k Δ[i]. Importantly, we only need the middle product here.
    prod = MP(p_tilde, S) .* Δ
    
    return prod
end
