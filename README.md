# Simple Polynomial Ring

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliahub.com/docs/SimplePolynomialRing/uWGaY)

Polynomial ring realization.

## Installation

```julia
(v1.8) pkg> add SimplePolynomialRing
```

This package supports Julia 1.8 and later.

## Usage

```julia
julia> using SimplePolynomialRing
```

## Construction

Polynomial could be constructed from Julia Expr or String.
```julia
julia> a = Pℤ{3, :x}(:(1 + x + 2x^2))
1 + x + 2x^2

julia> b = Pℤ{3, :x}("1 + x + 2x^2")
1 + x + 2x^2
```

As polynom's coefficients are modulo `P`, then the following polynomials are the same.
```julia
julia> c = Pℤ{3}(2)
2

julia> d = Pℤ{3}(5)
2

julia> c == d
true
```

```julia
julia> e = Pℤ₂([1,0,0,1,1,0,1]) # alias for `Pℤ{2, :x}`
1 + x^3 + x^4 + x^6
```

## Arithmetic

Operaions with two polynomials.
```julia
julia> p = Pℤ{3, :x}(:(1 + 2x + x^3 + 2x^5))
1 + 2x + x^3 + 2x^5

julia> q = Pℤ{3, :x}(:(1 + x + 2x^2))
1 + x + 2x^2

julia> p + q
2 + 2x^2 + x^3 + 2x^5

julia> p - q
x + x^2 + x^3 + 2x^5

julia> p * q
1 + x^2 + 2x^3 + x^4 + x^5 + 2x^6 + x^7

julia> p ÷ q
(2 + x + x^2 + x^3, 2 + 2x)

julia> p / q
2 + x + x^2 + x^3

julia> p % q
2 + 2x
```

Operations with a polynomial and a number.
```julia
julia> p = Pℤ{3, :x}(:(1 + 2x + x^3 + 2x^5))
1 + 2x + x^3 + 2x^5

julia> p * 5
2 + x + 2x^3 + x^5

julia> p ^ 5
1 + x + x^2 + x^3 + x^6 + 2x^7 + x^9 + 2x^10 + 2x^12 + x^14 + 2x^15 + 2x^16 + x^18 + 2x^20 + 2x^23 + 2x^25
```

## Some useful functions on polynomials

Extended GCD.
```julia
julia> p = Pℤ{3, :x}(:(1 + 2x + x^3 + 2x^5))
1 + 2x + x^3 + 2x^5

julia> q = Pℤ{3, :x}(:(1 + x + 2x^2))
1 + x + 2x^2

julia> gcdx(p, q)
(2, 2 + 2x, 2x^2 + 2x^3 + x^4)
```

Factorization and primality test.
```julia
julia> p = Pℤ{3, :x}(:(1 + 2x + x^3 + 2x^5))
1 + 2x + x^3 + 2x^5

julia> isprime(p)
false

julia> factor(p)
Dict{Pℤ{3, :x}, Int64} with 2 entries:
  2 + x              => 2
  1 + x + x^2 + 2x^3 => 1
```

Polynomial generation: all possible and prime.
```julia
julia> generate(3, 3)
54-element Vector{Pℤ{3, :x}}:
 x^3
 1 + x^3
 2 + x^3
 x + x^3
 1 + x + x^3
 2 + x + x^3
 2x + x^3
 1 + 2x + x^3
 2 + 2x + x^3
 x^2 + x^3
 1 + x^2 + x^3
 ⋮
 2 + 2x + x^2 + 2x^3
 2x^2 + 2x^3
 1 + 2x^2 + 2x^3
 2 + 2x^2 + 2x^3
 x + 2x^2 + 2x^3
 1 + x + 2x^2 + 2x^3
 2 + x + 2x^2 + 2x^3
 2x + 2x^2 + 2x^3
 1 + 2x + 2x^2 + 2x^3
 2 + 2x + 2x^2 + 2x^3

julia> primes(3, 3)
28-element Vector{Pℤ{3, :x}}:
 x
 1 + x
 2 + x
 2x
 1 + 2x
 2 + 2x
 1 + x^2
 2 + x + x^2
 2 + 2x + x^2
 2 + 2x^2
 1 + x + 2x^2
 ⋮
 1 + x + 2x^2 + x^3
 2 + 2x + 2x^2 + x^3
 1 + x + 2x^3
 2 + x + 2x^3
 2 + x^2 + 2x^3
 1 + x + x^2 + 2x^3
 2 + 2x + x^2 + 2x^3
 1 + 2x^2 + 2x^3
 2 + x + 2x^2 + 2x^3
 1 + 2x + 2x^2 + 2x^3

julia> irreducible(3, 3)
16-element Vector{Pℤ{3, :x}}:
 1 + 2x + x^3
 2 + 2x + x^3
 2 + x^2 + x^3
 2 + x + x^2 + x^3
 1 + 2x + x^2 + x^3
 1 + 2x^2 + x^3
 1 + x + 2x^2 + x^3
 2 + 2x + 2x^2 + x^3
 1 + x + 2x^3
 2 + x + 2x^3
 2 + x^2 + 2x^3
 1 + x + x^2 + 2x^3
 2 + 2x + x^2 + 2x^3
 1 + 2x^2 + 2x^3
 2 + x + 2x^2 + 2x^3
 1 + 2x + 2x^2 + 2x^3
```