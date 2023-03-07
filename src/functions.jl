# multithread=false
debug = false

multithread = true

include("types.jl")
using DSP

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
function MP(y::Array{<: Union{Integer, Rational, IntModQ}}, z::Array{<: Union{Integer, Rational, IntModQ}})
    # corollary 2 from https://hal.inria.fr/inria-00071921/document
    n = length(y)
    


    if n == 1
        return [y[1] * z[1]]
    elseif n == 2
        α = (z[1] + z[2]) * y[2]
        β = (y[1] - y[2]) * z[2]
        γ = (z[2] + z[3]) * y[1]
        return [α + β, γ - β]
    end


    # this is a work in progress
    if n > 12
        if typeof(y[1]) == IntModQ
            return IntModQ.(DSP.conv(convert.(Int,y),convert.(Int,z))[n:2*n-1])
        end
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



# performs the same computation as lagrange, but assumes G and H are precomputed
# use theorem 3.1 from https://dl.acm.org/doi/pdf/10.1145/120694.120697
function lagrange_precomputed(p::Vector{T1},G::Vector{T1},H::Vector{T1})::Vector{T1} where {T1 <: Union{Integer, Rational, IntModQ}}
    
    # TODO: the first half of the entries of G are zero... maybe make use of that?
    # also, we need straight up polynomial multiplication implemented
    n = length(p)
    
    
    # if typeof(p[1]) == IntModQ && n > 4
    #     d = n-1
    #     Gprime = G[1:d+2]
    #     Hprime = H[d+2:end]
    #     # Gprime =  T1[ (2(i & 1) - 1) * binomial(big(d + 1), big(i - 1)) for i = 1:d+2]
    #     # Hprime = T1[-(2((i-d-1) & 1) - 1) * binomial(big(-d - 1), big(i-d-1 - 1)) for i = d+2:(2d+4)]

    #     temp = DSP.conv(convert.(Int,Gprime),convert.(Int,p))[n:2*n]
    #     temp = convert.(Int,IntModQ.(temp))
    #     temp[1] = 0

    #     temp = IntModQ.(DSP.conv(convert.(Int,Hprime),temp))[2:n+1]
    #     return temp
    # end

    F1 = [0; MP([p; 0], G)[1:(end - 1)]]
    return MP(F1, H)[2:end]
end

@doc raw"""
    lagrange(p)

Lagrange interpolate a vector $p$ of polynomial outputs to a new position $degree+1$.

For a polynomial $P$ of degree at most d, apply Lagrange 
interpolation to the values $p = [P(0), P(1), ..., P(d)]$ 
to produce a vector of shifted values $[P(d+1), P(d + 2), ..., P(2d + 1)]$.

# Example
```jldoctest
julia> Main.OurModuleName.lagrange([1, -7, -31])
3-element Vector{Int64}:
  -71
 -127
 -199
```
"""
function lagrange(p::Vector{T1})::Vector{T1} where {T1 <: Union{Integer, Rational, IntModQ}}
    d = length(p) - 1
    G = T1[i < d+3 ? (2(i & 1) - 1) * binomial(big(d + 1), big(i - 1)) : 0 for i = 1:(2d+4)]
    H = T1[i > d+1 ? -(2((i-d-1) & 1) - 1) * binomial(big(-d - 1), big(i-d-1 - 1)) : 0 for i = 1:(2d+4)]
    return lagrange_precomputed(p,G,H)
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
function lagrange_matrix(matrices_to_interpolate::Vector{Array{T1, 2}}) where {T1 <: Union{Integer, Rational, IntModQ}}
    d = length(matrices_to_interpolate) - 1
    
    MATRIX_DIMENSION = size(matrices_to_interpolate[1], 1)

    # ith entry is a (matdim) by (matdim) matrix
    # I'd prefer to not initialize this to zeros but I'm not sure how else to do this
    # Todo: this should probably only be 3d + 1
    interpolated_matrices = [zeros(T1, MATRIX_DIMENSION, MATRIX_DIMENSION) for i = 1:3(d + 1)]

    
    G = T1[i < d+3 ? (2(i & 1) - 1) * binomial(big(d + 1), big(i - 1)) : 0 for i = 1:(2d+4)]
    H = T1[i > d+1 ? -(2((i-d-1) & 1) - 1) * binomial(big(-d - 1), big(i-d-1 - 1)) : 0 for i = 1:(2d+4)]


    if false#multithread
        Threads.@threads for thread_i = 0:MATRIX_DIMENSION*MATRIX_DIMENSION-1
            i = 1 + (thread_i ÷ MATRIX_DIMENSION)
            j = 1 + (thread_i - MATRIX_DIMENSION*(i-1))

            v1 = lagrange_precomputed([matrices_to_interpolate[k][i, j] for k in 1:d+1],G,H)
            v2 = lagrange_precomputed(v1, G, H)
            v3 = lagrange_precomputed(v2, G, H)
         
            for k in 1:d+1
                interpolated_matrices[k][i, j] = v1[k]
                interpolated_matrices[k + d + 1][i, j] = v2[k]
                interpolated_matrices[k + 2(d+1)][i, j] = v3[k]
            end
        end
    else
        for i in 1:MATRIX_DIMENSION
            for j in 1:MATRIX_DIMENSION
                v1 = lagrange_precomputed([matrices_to_interpolate[k][i, j] for k in 1:d+1],G,H)
                v2 = lagrange_precomputed(v1, G, H)
                v3 = lagrange_precomputed(v2, G, H)
             
                for k in 1:d+1
                    interpolated_matrices[k][i, j] = v1[k]
                    interpolated_matrices[k + d + 1][i, j] = v2[k]
                    interpolated_matrices[k + 2(d+1)][i, j] = v3[k]
                end
            end
        end
    end
    
    return interpolated_matrices
end

# need to nail down what the recursive product does

# IN: A list of polynomial matrices (M_ij evaluated at 0, 1, ..., max_degree); an integer a
# OUT: prod_{k=0}^{a maybe minus 1} M(k)

@doc raw"""
    matrix_product(T, a)

Compute the product T(a-1)T(a-2)...T(a) a matrix T.

# Example
```jldoctest
julia> Main.OurModuleName.matrix_product(BigInt.([0:2 3:5 6:8 ;;; 9:11 12:14 15:17 ]), 2^2)
3×3 Matrix{BigInt}:
  688662   762210   835758
  907092  1003968  1100844
 1125522  1245726  1365930
```
"""
function matrix_product(starter_matrices::Vector{Array{T, 2}}, a::S)::Array{T, 2} where {T <: Union{Integer, Rational, IntModQ}, S <: Integer}
    # currently assumes a is a power of 2 which is greater than d
    d = length(starter_matrices) - 1

    smat = copy(starter_matrices)
    num_steps = ceil(Int, log2(sqrt(a)))
    
    if debug
        println(num_steps)
        for j = 0:(num_steps - 1)
            println(j)
            smat = matrix_product_step(smat, d, j)
        end
    else
        for j = 0:(num_steps - 1)
            smat = matrix_product_step(smat, d, j)
        end
    end


    # then one big ol product

    #Multiply matrices according to price is right rules
    product_size = floor(Int, a / (2^num_steps))

    return prod(smat[1:product_size])

end



function matrix_product_step(smat::Vector{Array{T, 2}}, d::S, j::R)::Vector{Array{T, 2}} where {T <: Union{Integer, Rational, IntModQ}, S <: Integer, R <: Integer}
    
    smat = [smat ; lagrange_matrix(smat)]    
    upper_bound = (2^(j+1))*(d) + 1

    if multithread
        # quick and dirty way to initialize
        ans = smat[1:upper_bound]
        Threads.@threads for thread_i = 1:upper_bound
            ans[thread_i] = smat[2thread_i - 1] * smat[2thread_i]
        end
        return ans
    end

    return [smat[2i - 1] * smat[2i] for i = 1:upper_bound]
end
