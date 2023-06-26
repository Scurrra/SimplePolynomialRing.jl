# Combinatorics
"""
	Permutations{T}(a, t)

Returns permutations except that has first element of `a` at the end
"""
struct Permutations{T}
	a::T
	t::Int
end

Base.eltype(::Type{Permutations{T}}) where {T} = Vector{eltype(T)}

Base.length(p::Permutations) = length(p.a)^(p.t) - length(p.a)^(p.t-1)

function Base.iterate(p::Permutations, s=[ones(Int, p.t-1); 2])
	(!isempty(s) && s[end] > length(p.a) || p.t < 0) && return

	perm = p.a[s]
	s[1] += 1
	i = 1
	while i < length(s) && s[i] > length(p.a)
		s[i] = 1
		s[i+1] += 1
		i += 1
	end

	(perm, s)
end