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
function MP(y::Array{<: Union{Integer, Rational}}, z::Array{<: Union{Integer, Rational}})
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

Compute $\dfrac{1}{δ}$ where δ is the function from lemma 2 in [BGS07](https://specfun.inria.fr/bostan/publications/BoGaSc07.pdf).

For a degree $d$, return a vector $\left[\dfrac{1}{δ_0},..., \dfrac{1}{δ_d}\right]$,
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
function lagrange(p::Vector{T1}, a::T2; S::Vector=Vector(), δ::Vector=Vector(), Δ::Vector=Vector())::Vector{T1} where {T1 <: Union{Integer, Rational}, T2 <: Integer}
    d = length(p) - 1

    if a == d+1
        # use theorem 3.1 from https://dl.acm.org/doi/pdf/10.1145/120694.120697
        G = zeros(T1, 2d + 4)

        for i = 1:(d + 2)
            G[i] = (2(i & 1) - 1) * binomial(d + 1, i - 1)
            # G[i] = (-1)^(i-1) * binomial(d+1,i-1)
        end

        H = zeros(T1, 2d + 4)

        for i = 1:(d + 3) #I think 1:d+1 works too but I'm a little nervous about changing it
            H[d + 1 + i] = -(2(i & 1) - 1) * binomial(-d - 1, i - 1)
        end

        F1 = [0; MP([p; 0], G)[1:(end - 1)]]

        return MP(F1, H)[2:end]
    end

    if (length(δ) == 0)
        δ = compute_δ(d)
    end

    if (length(Δ) == 0)
        Δ = compute_Δ(a, d)
    end

    p_tilde = p .* δ

    if (length(S) == 0)
        S = zeros(Rational{Int}, 2d + 1)

        for i = 0:2d
            S[i + 1] = 1 // (a + i - d)
        end
    end

    #Q_k is the coefficient of x^{k + d} in ̃pS, and P(a + k) = Q_k Δ[i]. Importantly, we only need the middle product here.
    prod = MP(p_tilde, S) .* Δ
    
    return prod
end




# IN: list of d+1 matrices of polynomials of degree at most d
#     ie as a runs through the list, a_ij = P(0), P(1), ..., P(d)
# OUT: 3(d + 1) (?) matrices corresponding to P(d+1), ..., P(4d+1)

@doc raw"""
    lagrange_matrix(A)

Lagrange interpolate all entries of a list of $d+1$ matrices and return $3(d+1)$ matrices

Apply the following to each entry of the matrices, with $a = d+1, 2d+1, 3d+1$:
For a polynomial $P$, apply Lagrange interpolation to the values
$p = [P(0), P(1), ..., P(d)]$ to produce a vector of
shifted values $[P(a), P(a + 1), ..., P(a + d)]$.

# Example
```jldoctest
julia> Main.OurModuleName.lagrange_matrix([0:2 3:5 6:8 ;;; 9:11 12:14 15:17])
3×3×6 Array{Int64, 3}:
[:, :, 1] =
 18  21  24
 19  22  25
 20  23  26

[:, :, 2] =
 27  30  33
 28  31  34
 29  32  35

[:, :, 3] =
 36  39  42
 37  40  43
 38  41  44

[:, :, 4] =
 45  48  51
 46  49  52
 47  50  53

[:, :, 5] =
 54  57  60
 55  58  61
 56  59  62

[:, :, 6] =
 63  66  69
 64  67  70
 65  68  71
```
"""
function lagrange_matrix(matrices_to_interpolate::Array{T, 3}) where {T <: Union{Integer, Rational}}
    d = size(matrices_to_interpolate, 3) - 1

    MATRIX_DIMENSION = size(matrices_to_interpolate, 1)

    # ith entry is a (matdim) by (matdim) matrix
    # I'd prefer to not initialize this to zeros but I'm not sure how else to do this
    interpolated_matrices = zeros(T, MATRIX_DIMENSION, MATRIX_DIMENSION, 3(d+1))
    

    δ = compute_δ(d)
    Δ2 = compute_Δ(2d+2, d)
    Δ3 = compute_Δ(3d+3, d)

    S2 = zeros(Rational{Int}, 2d + 1)
    S3 = zeros(Rational{Int}, 2d + 1)
    for i = 0:2d
        S2[i + 1] = 1 // (2d+2 + i - d)
        S3[i + 1] = 1 // (3d+3 + i - d)
    end

    # the thing to do would be to compute \delta here and pass it in to each call
    for i in 1:MATRIX_DIMENSION
        for j in 1:MATRIX_DIMENSION
            v1 = lagrange([matrices_to_interpolate[i, j, k] for k in 1:d+1], d+1)
            v2 = lagrange([matrices_to_interpolate[i, j, k] for k in 1:d+1], 2d+2; S = S2, δ, Δ = Δ2)
            v3 = lagrange([matrices_to_interpolate[i, j, k] for k in 1:d+1], 3d+3; S = S3, δ, Δ = Δ3)

            for k in 1:d+1
                interpolated_matrices[i, j, k] = v1[k]
                interpolated_matrices[i, j, k + d + 1] = v2[k]
                interpolated_matrices[i, j, k + 2(d+1)] = v3[k]
            end
        end
    end
    
    return interpolated_matrices
end

# need to nail down what the recursive product does

# IN: A list of polynomial matrices (M_ij evaluated at 0, 1, ..., max_degree); an integer a
# OUT: prod_{k=0}^{a maybe minus 1} M(k)
function matrix_product(starter_matrices::Array{T, 3}, a::S) where {T <: Union{Integer, Rational}, S<:Integer}
    d = size(starter_matrices, 3) - 1
    # first, need to lagrange interpolate to get the nearest power of 2 matrices
    # next, do the recursive product formula
end