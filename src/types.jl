import Base: +, -, *, ^, zero, convert, promote_rule, binomial

characteristic_q = 3

struct IntModQ <: Number
    n::Int

    IntModQ(n) = new(mod(n, characteristic_q))
end

+(a::IntModQ, b::IntModQ) = IntModQ(a.n + b.n)
-(a::IntModQ, b::IntModQ) = IntModQ(a.n - b.n)
*(a::IntModQ, b::IntModQ) = IntModQ(a.n * b.n)
^(a::IntModQ, b::Int) = IntModQ(a.n ^ b)

function square(a::IntModQ)
    return IntModQ(a.n * a.n)
end

zero(::Type{IntModQ}) = IntModQ(0)

convert(::Type{IntModQ}, x::T) where {T <: Integer} = IntModQ(x)
convert(::Type{T}, x::IntModQ) where {T <: Integer} = T(x.n)
promote_rule(::Type{IntModQ}, ::Type{<: Integer}) = IntModQ

# function binomial(n::IntModQ, k::IntModQ)::IntModQ
#     #Todo: improve this mod p
#     return binomial(convert(Int, n), convert(Int, k))
# end