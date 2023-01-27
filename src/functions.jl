@doc raw"""
    MP(y, z)

Compute the middle product of polynomials $y$ and $z$.

For vectors $y = [y_1,...,y_n]$ and $z = [z_1,...,z_{2n-1}]$ of polynomial coefficients,
return $[c_{n-1},..., c_{2n-1}]$, where $c_k = \sum_{i + j = k} a_i b_j$.

# Examples
```jldoctest
julia> Main.OurModuleName.MP([4, 6, 38, -1], [1, 2, 3, 4, 8, 9, 26, 7])
4-element Vector{Int64}:
 109
 168
 233
 458
```

```jldoctest
julia> Main.OurModuleName.MP([-7//10, 1, -2], [9//4, -1//4, -2//5, 9//8, -4//9, -1//2])
3-element Vector{Rational{Int64}}:
 -447//100
  -11//16
  161//72
```
"""
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



@doc raw"""
    compute_δ(d)

Compute the δ function from lemma 2 in [BGS07](https://specfun.inria.fr/bostan/publications/BoGaSc07.pdf).

For a degree $d$, return a vector $[δ_0,..., δ_d]$,
where $δ_i = \prod_{j = 0, j \neq i}^d (i - j)$.

# Example
```jldoctest
julia> Main.OurModuleName.compute_δ(5)
6-element Vector{Rational}:
 -1//120
  1//24
 -1//12
  1//12
 -1//24
  1//120
```
"""
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



@doc raw"""
    compute_Δ(a, d)

Compute the Δ function from lemma 2 in [BGS07](https://specfun.inria.fr/bostan/publications/BoGaSc07.pdf).

For a degree $d$ and shift $a > d$, return a vector $[Δ_0,..., Δ_d]$,
where $Δ_i = \prod_{j = 0}^d (a + i - j)$.

# Example
```jldoctest
julia> Main.OurModuleName.compute_Δ(5, 3)
4-element Vector{Integer}:
  120
  360
  840
 1680
```
"""
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



@doc raw"""
    lagrange(p, a)

Lagrange interpolate a vector $p$ of polynomial outputs to a new position $a$.

For a polynomial $P$, apply Lagrange interpolation to the values
$p = [P(0), P(1), ..., P(d)]$ to produce a vector of
shifted values $[P(a), P(a + 1), ..., P(a + d)]$.

# Example
```jldoctest
julia> Main.OurModuleName.lagrange([1//4, -7//4, -31//4], 3)
3-element Vector{Rational{Int64}}:
  -71//4
 -127//4
 -199//4
```
"""
function lagrange(p::Vector{<: Union{Integer, Rational}}, a::T) where {T <: Integer}
    d = length(p) - 1

    if a == d+1
        # use theorem 3.1 from https://dl.acm.org/doi/pdf/10.1145/120694.120697
        G = fill(0//1, 2d + 4)

        for i = 1:(d + 2)
            G[i] = (2(i & 1) - 1) * binomial(d + 1, i - 1)
            # G[i] = (-1)^(i-1) * binomial(d+1,i-1)
        end

        H = fill(0//1,2d + 4)

        for i = 1:(d + 3) #I think 1:d+1 works too but I'm a little nervous about changing it
            H[d + 1 + i] = -(2(i & 1) - 1) * binomial(-d - 1, i - 1)
        end

        F1 = [0; MP([p; 0], G)[1:(end - 1)]]

        return MP(F1, H)[2:end]
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