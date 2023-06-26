module SimplePolynomialRing

using SparseArrays

import Primes: isprime, factor, primes
import Base: (+), (-), (*), (^), (÷), (%), (/)

include("tools.jl")

"""
	Pℤ{P, X}

Polynomial of ℤₚ[X].

# Examples
```jldoctest
julia> a = Pℤ{3, :x}(:(1 + x + 2x^2))
1 + x + 2x^2

julia> b = Pℤ{3, :x}("1 + x + 2x^2")
1 + x + 2x^2

julia> c = Pℤ{3}(2)
2

julia> d = Pℤ{3}(5)
2

julia> e = Pℤ₂([1,0,0,1,1,0,1])
1 + x^3 + x^4 + x^6
```
"""
struct Pℤ{P,X}
    coeffs::SparseVector{Int,Int}

    function Pℤ{P,X}(coeffs::SparseVector{Int,Int}) where {P,X}
        (typeof(P) != Int || P < 2) && throw("`P` must be natural number")
        typeof(X) != Symbol && throw("`X` must be Symbol")

        for d in coeffs.nzind
            if coeffs[d] < 0
                while coeffs[d] < 0
                    coeffs[d] += P
                end
            else
                coeffs[d] %= P
            end
        end

        dropzeros!(coeffs)
        # coeffs = sparsevec(coeffs.nzind, coeffs.nzval)
        # coeffs = length(coeffs.nzind)==0 ? sparsevec([0]) : coeffs
        new{P,X}(coeffs)
    end

    function Pℤ{P,X}(coeffs::Vector{Int}) where {P,X}
        coeffs = sparsevec(coeffs)
        Pℤ{P,X}(coeffs)
    end

    function Pℤ{P,X}(poly::Int) where {P,X}
        Pℤ{P,X}(sparsevec([(poly + P) % P]))
    end

    function Pℤ{P,X}(poly::Symbol) where {P,X}
        Pℤ{P,poly}(sparsevec([0, 1]))
    end

    function Pℤ{P,X}(poly::Expr) where {P,X}
        (typeof(P) != Int || P < 2) && throw("`P` must be natural number")
        typeof(X) != Symbol && throw("`X` must be Symbol")

        coeffs = Dict{Int,Int}()
        if length(poly.args) == 0
            coeffs[1] = 0
        elseif poly.args[1] == :+
            for arg in poly.args[2:end]
                if typeof(arg) == Int
                    coeffs[1] = arg
                elseif typeof(arg) == Symbol
                    coeffs[2] = 1
                elseif typeof(arg) == Expr
                    inner = arg.args[2:end]
                    if typeof.(inner) == [Symbol, Int]
                        coeffs[inner[2]+1] = 1
                    elseif typeof.(inner) == [Int, Expr]
                        coeffs[inner[2].args[3]+1] = inner[1]
                    elseif typeof.(inner) == [Int, Symbol]
                        coeffs[2] = inner[1]
                    end
                end
            end
        elseif poly.args[1] == :^
            coeffs[poly.args[end]+1] = 1
        elseif poly.args[1] == :*
            if typeof(poly.args[end]) == Expr
                coeffs[poly.args[end].args[3]+1] = poly.args[2]
            elseif typeof(poly.args[end]) == Symbol
                coeffs[2] = poly.args[2]
            end
        end

        Pℤ{P,X}(coeffs |> sparsevec)
    end

    function Pℤ{P,X}(poly::String) where {P,X}
        Pℤ{P,X}(poly |> Meta.parse)
    end

    Pℤ{P}(poly::Union{Expr,String,Int,Symbol,Vector{Int}}) where {P} = Pℤ{P,:x}(poly)
end
Base.one(::Pℤ{P,X}) where {P,X} = Pℤ{P,X}(1)
Base.zero(::Pℤ{P,X}) where {P,X} = Pℤ{P,X}(0)
Base.one(::Type{Pℤ{P,X}}) where {P,X} = Pℤ{P,X}(1)
Base.zero(::Type{Pℤ{P,X}}) where {P,X} = Pℤ{P,X}(0)
function Base.getproperty(p::Pℤ{P,X}, name::Symbol) where {P,X}
    if name == :P
        return P
    elseif name == :X
        return X
    else
        return getfield(p, name)
    end
end

Base.copy(p::Pℤ{P,X}) where {P,X} = Pℤ{P,X}(p.coeffs)
Base.hash(p::Pℤ{P,X}, h::UInt) where {P,X} = hash((p.coeffs, p.P, p.X), h)

"""
    degree(p::Pℤ)

Return the degree of the polynomial `p`.
"""
degree(p::Pℤ{P,X}) where {P,X} = length(p.coeffs.nzind) == 0 ? 0 : (p.coeffs.nzind |> maximum) - 1

import Base: (==), (!=)
Base.isequal(left::Pℤ{P,X}, right::Pℤ{P,X}) where {P,X} = left.coeffs.nzind == right.coeffs.nzind && left.coeffs.nzval == right.coeffs.nzval
(==)(left::Pℤ{P,X}, right::Pℤ{P,X}) where {P,X} = isequal(left, right)
(==)(v::Vector{Pℤ{P,X}}, p::Pℤ{P,X}) where {P,X} = broadcast(==, v, p)
(==)(p::Pℤ{P,X}, v::Vector{Pℤ{P,X}}) where {P,X} = broadcast(==, p, v)

Base.Vector(p::Pℤ{P,X}) where {P,X} = length(p.coeffs.nzind) == 0 ? [0] :
                                      let v = zeros(Int, maximum(p.coeffs.nzind[p.coeffs.nzval.!=0]))
    v[p.coeffs.nzind[p.coeffs.nzval.!=0]] .= p.coeffs.nzval[p.coeffs.nzval.!=0]
    return v
end

function Base.unique(v::Vector{Pℤ{P,X}}) where {P,X}
    uv = []
    while length(v) > 0
        push!(uv, v[1])
        v = v[.!(v == v[1])]
    end
    sort!(uv, by=x -> degree(x))
end

Base.Symbol(p::Pℤ) = Symbol(string(p))

Base.size(::Pℤ) = ()
Base.getindex(p::Pℤ, i) = p
Base.Broadcast.broadcastable(p::Pℤ) = Ref(p)

function Base.show(io::IO, p::Pℤ{P,X}) where {P,X}
    degrees = p.coeffs.nzind
    if length(degrees) == 0
        return print(io, "0")
    end
    p_str = ""

    for i in 1:length(degrees)-1
        if degrees[i] == 1 # x ^ 0
            p_str *= "$(p.coeffs[1])"
        elseif degrees[i] == 2 # x ^ 1
            if p.coeffs[2] != 1
                p_str *= "$(p.coeffs[2])"
            end
            p_str *= String(X)
        else
            if p.coeffs[degrees[i]] != 1
                p_str *= "$(p.coeffs[degrees[i]])"
            end
            p_str *= String(X) * "^$(degrees[i]-1)"
        end
        p_str *= " + "
    end

    if degrees[end] == 1
        p_str *= "$(p.coeffs[1])"
    elseif degrees[end] == 2 # x ^ 1
        if p.coeffs[2] != 1
            p_str *= "$(p.coeffs[2])"
        end
        p_str *= String(X)
    else
        if p.coeffs[degrees[end]] != 1
            p_str *= "$(p.coeffs[degrees[end]])"
        end
        p_str *= String(X) * "^$(degrees[end]-1)"
    end

    print(io, p_str)
end

function +(left::Pℤ{P,X}, right::Pℤ{P,X})::Pℤ{P,X} where {P,X}
    left_length = left.coeffs |> length
    right_length = right.coeffs |> length
    ans_length = max(left_length, right_length)

    Pℤ{P,X}(
        [left.coeffs; spzeros(Int, ans_length - left_length)] +
        [right.coeffs; spzeros(Int, ans_length - right_length)]
    )
end

function -(left::Pℤ{P,X}, right::Pℤ{P,X})::Pℤ{P,X} where {P,X}
    left_length = left.coeffs |> length
    right_length = right.coeffs |> length
    ans_length = max(left_length, right_length)

    Pℤ{P,X}(
        [left.coeffs; spzeros(Int, ans_length - left_length)] -
        [right.coeffs; spzeros(Int, ans_length - right_length)]
    )
end

(*)(left::Pℤ{P,X}, right::Int) where {P,X} = degree(left) == 0 ? left : Pℤ{P,X}(left.coeffs .* right)

(*)(left::Int, right::Pℤ{P,X}) where {P,X} = right * left

function *(left::Pℤ{P,X}, right::Pℤ{P,X}) where {P,X}
    left_degrees = left.coeffs.nzind
    right_degrees = right.coeffs.nzind

    coeffs = spzeros(Int, degree(left) + degree(right) + 1)

    for ld in left_degrees
        for rd in right_degrees
            coeffs[ld+rd-1] +=
                left.coeffs[ld] * right.coeffs[rd]
        end
    end
    Pℤ{P,X}(coeffs)
end

(^)(left::Pℤ{P,X}, right::Int) where {P,X} = Base.power_by_squaring(left, right)

function ÷(left::Pℤ{P,X}, right::Pℤ{P,X}) where {P,X}
    div = Pℤ{P,X}(:())
    rc = right.coeffs[degree(right)+1]
    while degree(left) >= degree(right) && left != Pℤ{P,X}(:())
        c = 1
        while c * rc % P != left.coeffs[degree(left)+1]
            c += 1
        end

        div += Pℤ{P,X}("$(c)$(X)^$(degree(left)-degree(right))")
        left -= right * Pℤ{P,X}("$(c)$(X)^$(degree(left)-degree(right))")
    end
    (div, left)
end

function /(left::Pℤ{P,X}, right::Pℤ{P,X}) where {P,X}
    return ÷(left, right)[1]
end
function %(left::Pℤ{P,X}, right::Pℤ{P,X}) where {P,X}
    return ÷(left, right)[2]
end

+(v::Vector{Pℤ{P,X}}, p::Pℤ{P,X}) where {P,X} = broadcast(+, v, p)
+(p::Pℤ{P,X}, v::Vector{Pℤ{P,X}}) where {P,X} = broadcast(+, p, v)
-(v::Vector{Pℤ{P,X}}, p::Pℤ{P,X}) where {P,X} = broadcast(-, v, p)
-(p::Pℤ{P,X}, v::Vector{Pℤ{P,X}}) where {P,X} = broadcast(-, p, v)
*(v::Vector{Pℤ{P,X}}, p::Pℤ{P,X}) where {P,X} = broadcast(*, v, p)
*(p::Pℤ{P,X}, v::Vector{Pℤ{P,X}}) where {P,X} = broadcast(*, p, v)
*(v::Vector{Pℤ{P,X}}, n::Int) where {P,X} = broadcast(*, v, n)
^(v::Vector{Pℤ{P,X}}, n::Int) where {P,X} = broadcast(^, v, n)
%(v::Vector{Pℤ{P,X}}, p::Pℤ{P,X}) where {P,X} = broadcast(%, v, p)
%(p::Pℤ{P,X}, v::Vector{Pℤ{P,X}}) where {P,X} = broadcast(%, p, v)
/(v::Vector{Pℤ{P,X}}, p::Pℤ{P,X}) where {P,X} = broadcast(/, v, p)
/(p::Pℤ{P,X}, v::Vector{Pℤ{P,X}}) where {P,X} = broadcast(/, p, v)

"""
    gcdx(a::Pℤ{P,X}, b::Pℤ{P,X})

Computes the greatest common divisor of `a` and `b` and their Bézout coefficients, i.e. the polynomials `u` and `v` that satisfy `ua+vb = d = gcd(a, b)`. `gcdx(a, b) returns `(d, u, v)`.
"""
function Base.gcdx(a::Pℤ{P,X}, b::Pℤ{P,X}) where {P,X}
    if a == Pℤ{P,X}(0)
        return b, Pℤ{P,X}(0), Pℤ{P,X}(1)
    end

    d, x₁, y₁ = gcdx(b % a, a)
    return d, y₁ - (b / a) * x₁, x₁
end

# Factorization
struct PℤGenerator{P,X}
    a::UnitRange
    t::Int

    function PℤGenerator{P,X}(n::Int) where {P,X}
        (typeof(P) != Int || P < 2) && throw("`P` must be natural number")
        typeof(X) != Symbol && throw("`X` must be Symbol")
        n < 0 && throw("`n` (degree of polynomial) must be non negative number")

        new{P,X}(0:(P-1), n + 1)
    end
end

Base.eltype(::Type{PℤGenerator{P,X}}) where {P,X} = Pℤ{P,X}

Base.length(p::PℤGenerator) = p.t == 1 ? length(p.a) : length(p.a)^(p.t) - length(p.a)^(p.t - 1)

function Base.iterate(p::PℤGenerator{P,X}, s=[ones(Int, p.t - 1); [1, 2][Int(p.t != 1)+1]]) where {P,X}
    (!isempty(s) && s[end] > length(p.a) || p.t < 0) && return

    perm = p.a[s]
    s[1] += 1
    i = 1
    while i < length(s) && s[i] > length(p.a)
        s[i] = 1
        s[i+1] += 1
        i += 1
    end

    (Pℤ{P,X}(perm), s)
end

"""
    generate(n::Int, p::Int; x::Symbol=:x)

Generates all the polynomials over ℤₚ up to degree `n`.
"""
function generate(n::Int, p::Int; x::Symbol=:x)
    n < 0 && throw("`n` (degree of polynomial) must be non negative number")
    p < 1 && throw("`p` must be natural number")

    n == 0 ? Pℤ{p,x}.(0:(p-1)) : Pℤ{p,x}.(Permutations(0:(p-1), n + 1))
end

"""
    primes(n::Int, p::Int; x::Symbol=:x)

Returns a collection of the prime polynomials over ℤₚ up to degree `n`.
"""
function primes(n::Int, p::Int; x::Symbol=:x)
    !isprime(p) && throw("`p` must be prime number")
    ZERO = Pℤ{p,x}(0)
    primes_list = generate(1, p, x=x) |> collect

    for i in 2:n
        gen_list = generate(i, p, x=x) |> collect

        append!(primes_list, gen_list[sum([
            gen_list % primes_list[j] == ZERO
            for j in 1:round(Int, length(primes_list) / 2)
        ]).==0])
    end

    primes_list
end

"""
    factor(p::Pℤ)

Compute the prime factorization of a polynomial `p`.
"""
function factor(p::Pℤ{P,X}) where {P,X}
    (p == Pℤ{P,X}(1) || p == Pℤ{P,X}(0)) && return Dict(p => 1)

    !isprime(P) && throw("`P` must be prime number")
    ZERO = Pℤ{P,X}(0)

    rainbow = primes(round(Int, degree(p) / 2), P; x=X)

    factorization = Dict{Pℤ{P,X},Int}()

    for r in rainbow
        if p % r == ZERO
            factorization[r] = 0
            while p % r == ZERO
                factorization[r] += 1
                p = p / r
            end
        end
    end

    if p != Pℤ{P,X}(1)
        factorization[p] = 1
    end
    factorization
end

"""
    isprime(p::Pℤ) -> Bool

Returns `true` if polynomial `p` is prime, and `false` otherwise.
"""
isprime(p::Pℤ) = length(factor(p)) == 1

"""
    irreducible(n::Int, p::Int; x::Symbol=:x)

Generates irreducible polynomials over ℤₚ of degree `n`.
"""
function irreducible(n::Int, p::Int; x::Symbol=:x)
    n < 1 && throw("`n` (degree of polynomial) must be natural number")
    p < 1 && throw("`p` must be natural number")
    !isprime(p) && throw("`p` must be prime number")

    filter(
        p -> degree(p) == n,
        primes(n, p; x=x)
    )
end

"""
	Pℤ₂

Alias for Pℤ{2, :x}
"""
Pℤ₂ = Pℤ{2}

export Pℤ, Pℤ₂, degree, generate, primes, isprime, factor, irreducible

end # module PolynomialRing
