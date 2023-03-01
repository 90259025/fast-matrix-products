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
    # corollary 2 from https://hal.inria.fr/inria-00071921/document
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
6-element Vector{Rational{BigInt}}:
 -1//120
  1//24
 -1//12
  1//12
 -1//24
  1//120
```
"""
function compute_δ(d::T) where T <: Integer
    δ = Array{Rational{BigInt}}(undef, d + 1)
    
    if iseven(d)
        δ[1] = 1 // factorial(big(d)) 
    else
        δ[1] = -1 // factorial(big(d))
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
4-element Vector{BigInt}:
  120
  360
  840
 1680
```
"""
function compute_Δ(a::T, d::S) where {T <: Integer, S <: Integer}
    Δ = Array{BigInt}(undef, d + 1)
    Δ₀ = big(1)
    
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
            G[i] = (2(i & 1) - 1) * binomial(big(d + 1), big(i - 1))
            # G[i] = (-1)^(i-1) * binomial(d+1,i-1)
        end

        H = zeros(T1, 2d + 4)

        for i = 1:(d + 3) #I think 1:d+1 works too but I'm a little nervous about changing it
            H[d + 1 + i] = -(2(i & 1) - 1) * binomial(big(-d - 1), big(i - 1))
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
        S = zeros(Rational{BigInt}, 2d + 1)

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
    # Todo: this should probably only be 3d + 1
    interpolated_matrices = zeros(T, MATRIX_DIMENSION, MATRIX_DIMENSION, 3(d + 1))
    

    # δ = compute_δ(d)
    # Δ2 = compute_Δ(2d+2, d)
    # Δ3 = compute_Δ(3d+3, d)

    # S2 = zeros(Rational{BigInt}, 2d + 1)
    # S3 = zeros(Rational{BigInt}, 2d + 1)
    # for i = 0:2d
    #     S2[i + 1] = 1 // (2d+2 + i - d)
    #     S3[i + 1] = 1 // (3d+3 + i - d)
    # end

    # the thing to do would be to compute \delta here and pass it in to each call
    for i in 1:MATRIX_DIMENSION
        for j in 1:MATRIX_DIMENSION
            v1 = lagrange([matrices_to_interpolate[i, j, k] for k in 1:d+1], d+1)
            # v2 = lagrange(v1, d+1)
            # v3 = lagrange(v2, d+1)
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
function matrix_product(starter_matrices::Array{T, 3}, a::S) where {T <: Union{Integer, Rational}, S <: Integer}
    # currently assumes a is a power of 2 which is greater than d
    d = size(starter_matrices, 3) - 1

    smat = copy(starter_matrices)
    num_steps = ceil(Int, log2(sqrt(a)))

    for j = 0:(num_steps - 1)
        smat = matrix_product_step(smat, d, j)
    end

    # then one big ol product
    smat_prod = smat[:, :, 1]

    #Multiply matrices according to price is right rules
    product_size = floor(Int, a / (2^num_steps))

    # println(product_size, smat)

    for j = 2:product_size
        smat_prod = smat_prod * smat[:, :, j]
    end
    
    return smat_prod
end



function matrix_product_step(smat::Array{T, 3}, d::S, j::R) where {T <: Union{Integer, Rational}, S <: Integer, R <: Integer}
    smat = [smat ;;; lagrange_matrix(smat)]

    # println("before slicing: ", smat)

    # I'm sure this can be cleaned up
    for i in 1:(2^(j+1))*(d) + 1
        smat[:, :, i] = smat[:, :, 2i - 1] * smat[:, :, 2i]
    end

    smat = smat[:, :, 1:(2^(j+1))*(d) + 1]

    return smat
end

# println(mod.(matrix_product(BigInt.([9 9 -4; 6 1 1; 2 2 -8 ;;; 0 4 -20; -6 13 6; -1 -12 9 ;;; -25 3 -48; -22 43 11; 10 -40 42 ]), 2),2053)  == BigInt.([2003 201 1891; 2046 25 1948; 2049 130 1953]))

# println(matrix_product(BigInt.([-7 -8 -7; 0 6 3; -8 -6 -1 ;;; -13 -9 -16; -3 7 19; -6 -10 7 ;;; -27 -10 -23; 4 16 53; -12 -4 15 ]), 2)  == BigInt.([157 77 -89; -36 12 135; 128 40 7]))
# println(matrix_product(BigInt.([1 -3 0; 1 -1 8; -1 -2 4 ;;; 3 -12 9; 8 1 10; -4 -6 -7 ;;; 9 -31 24; 27 -3 14; -17 -22 -38 ]), 16)  == BigInt.([26232097201847637501656274121044930537 -252761813710026155381080448027519775048 -252158182194751717460963166378577530153; 87931715687714391870923057642403710511 -991443661719657974170019250074987959386 -1066541778933255595259716322568147718527; 49033727672125983633121590997052304384 -1118752398295368069898807765066599654123 -1223672762795933290451505791107355310744]))
# println(matrix_product(BigInt.([1 -3 0; 1 -1 8; -1 -2 4 ;;; 3 -12 9; 8 1 10; -4 -6 -7 ;;; 9 -31 24; 27 -3 14; -17 -22 -38 ]), 16))

# println(mod.(matrix_product(BigInt.([4 4 -5 4 1; -9 -10 9 -2 0; -3 -4 -5 5 -3; -9 4 -7 4 -1; -9 -8 0 -9 1 ;;; 3 13 1 1 -8; -8 0 13 -4 7; -22 -1 -17 5 7; -16 18 -22 1 -10; 0 -14 -13 -9 2 ;;; -18 124 31 -156 15; -11 140 39 26 106; -93 126 5 -81 127; -73 150 -79 34 -71; 97 -122 -118 33 -141 ;;; -185 541 103 -833 232; -78 734 99 136 519; -240 635 157 -433 591; -366 676 -238 235 -220; 438 -560 -531 171 -764 ;;; -720 1588 235 -2612 949; -341 2322 205 398 1612; -463 1928 583 -1327 1753; -1225 2064 -583 856 -493; 1251 -1700 -1612 483 -2395 ]), 256),2053)  == BigInt.([30 1832 46 1524 704; 1950 975 1437 1939 1120; 1116 1999 1033 811 773; 1635 1107 642 63 685; 961 525 1309 1212 2006]))


if test_δ
    println("Testing delta...")
    # Test delta
    println(compute_δ(5) == [-1//120, 1//24, -1//12, 1//12, -1//24, 1//120])
    println(compute_δ(30)==[1//265252859812191058636308480000000, -1//8841761993739701954543616000000, 1//609776689223427721003008000000, -1//65333216702510112964608000000, 1//9678995067038535254016000000, -1//1861345205199718318080000000, 1//446722849247932396339200000, -1//130294164363980282265600000, 1//45319709343993141657600000, -1//18539881095269921587200000, 1//8828514807271391232000000, -1//4855683143999265177600000, 1//3066747248841641164800000, -1//2214873013052296396800000, 1//1824013069572479385600000, -1//1710012252724199424000000, 1//1824013069572479385600000, -1//2214873013052296396800000, 1//3066747248841641164800000, -1//4855683143999265177600000, 1//8828514807271391232000000, -1//18539881095269921587200000, 1//45319709343993141657600000, -1//130294164363980282265600000, 1//446722849247932396339200000, -1//1861345205199718318080000000, 1//9678995067038535254016000000, -1//65333216702510112964608000000, 1//609776689223427721003008000000, -1//8841761993739701954543616000000, 1//265252859812191058636308480000000])
end

if test_Δ
    println("Testing Delta...")
    println(compute_Δ(5,3)==[120, 360, 840, 1680])
    # This overflows, TODO fix
    println(compute_Δ(50,30)==BigInt.([250023166568122405335142878324631535616000000000, 637559074748712133604614339727810415820800000000, 1578717708901572902259045031706959124889600000000, 3803274480535607446351335758203128800870400000000, 8929427041257513134911831780129085010739200000000, 20463270302881800934172947829462486482944000000000, 45837725478455234092547403137995969721794560000000, 100490398164305705510584691494837318236241920000000, 215868262723323367393107855803724609544519680000000, 454865267881288524149762981872133998683095040000000, 941100554237148670654682031459587583482265600000000, 1913571126948868963664520130634494753080606720000000, 3827142253897737927329040261268989506161213440000000, 7534686312361171544429048014373323090254888960000000, 14612724969427726631619971906663414478070087680000000, 27936091853317712678097005115680057090428108800000000, 52679487494827686764411495360996679084807290880000000, 98042379504262639255988060810743819407835791360000000, 180185994764590796470464544192718370803590103040000000, 327179832598862235696369830244672831195992555520000000, 587245853382573243557586874798130722659473817600000000, 1042361389754067507314716702766682032720566026240000000, 1830488294202264890894136648761002594045872046080000000, 3181562987542031834173142270465552127746396651520000000, 5475247932049078040204942511963973429144961679360000000, 9332808975083655750349333827211318345133457408000000000, 15762077380141285267256652685956893205114283622400000000, 26384346918932151425625266452580016886821735628800000000, 43786788503759740663803633687260453556853093171200000000, 72065756079104573175843480443616163145654049177600000000, 117658377272007466409540376234475368401067835392000000000]))
end