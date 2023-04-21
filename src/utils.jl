# Copyright (c) 2021, 2023 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

function _isisomorphic(A::fq_nmod_mat, B::fq_nmod_mat)
    n = ncols(A)
    nrA = nrows(A)
    nrB = nrows(B)
    F = base_ring(A)
    V = VectorSpace(F, n)
    AVS, _ = sub(V, [V(A[i, :]) for i in 1:nrA])
    BVS, _ = sub(V, [V(B[i, :]) for i in 1:nrB])
    return is_isomorphic(AVS, BVS)
end

"""
    reverse(v::fq_nmod_mat)
    reverse!(v::fq_nmod_mat)

Return the reverse of the vector `v`.
"""
function reverse(v::fq_nmod_mat)
    nr, nc = size(v)
    u = deepcopy(v)
    if nr == 1
        left = 1
        right = nc
        while left < right
            temp = u[1, left]
            u[1, left] = u[1, right]
            u[1, right] = temp
            left += 1
            right -= 1
        end
        return u
    elseif nc == 1
        left = 1
        right = nr
        while left < right
            temp = u[left, 1]
            u[left, 1] = u[right, 1]
            u[right, 1] = temp
            left += 1
            right -= 1
        end
        return u
    else
        throw(ArgumentError("Matrix must be a vector."))
    end
end

function reverse!(v::fq_nmod_mat)
    nr, nc = size(v)
    if nr == 1      
        left = 1
        right = nc
        while left < right
            temp = v[1, left]
            v[1, left] = v[1, right]
            v[1, right] = temp
            left += 1
            right -= 1
        end
    elseif nc == 1
        left = 1
        right = nr
        while left < right
            temp = v[left, 1]
            v[left, 1] = v[right, 1]
            v[right, 1] = temp
            left += 1
            right -= 1
        end
    else
        throw(ArgumentError("Matrix must be a vector."))
    end
end

"""
    circshift(v::fq_nmod_mat, l::Int)

Return the circular shift of the vector `v` by `l` bits.

This is an overload of Base.circshift for type `fq_nmod_mat`.
Either the number of rows or the number of columns must have dimension one.1
"""
function circshift(v::fq_nmod_mat, l::Int)
    nr, nc = size(v)
    if nr == 1
        l = l % nc
        l < 0 && (l = nc + l;)
        vshift = zero_matrix(base_ring(v), 1, nc)
        vshift[1, 1:l] = v[1, nc - l + 1:nc]
        vshift[1, l + 1:end] = v[1, 1:nc - l]
    elseif nc == 1
        l = l % nr
        l < 0 && (l = nr + l;)
        vshift = zero_matrix(base_ring(v), nr, 1)
        vshift[1:l, 1] = v[nr - l + 1:nr, 1]
        vshift[l + 1:end, 1] = v[1:nr - l, 1]
    else
        throw(ArgumentError("Input matrix must be a vector."))
    end
    return vshift
end

"""
    ⊕(A::fq_nmod_mat, B::fq_nmod_mat)
    directsum(A::fq_nmod_mat, B::fq_nmod_mat)

Return the direct sum of the two matrices `A` and `B`.
"""
function ⊕(A::fq_nmod_mat, B::fq_nmod_mat)
    base_ring(A) == base_ring(B) || error("Matrices must be over the same base ring in directsum.")

    return vcat(hcat(A, zero_matrix(base_ring(B), nrows(A), ncols(B))),
        hcat(zero_matrix(base_ring(A), nrows(B), ncols(A)), B))
end
directsum(A::fq_nmod_mat, B::fq_nmod_mat) = A ⊕ B

"""
    ⊗(A::fq_nmod_mat, B::fq_nmod_mat)
    kron(A::fq_nmod_mat, B::fq_nmod_mat)
    tensorproduct(A::fq_nmod_mat, B::fq_nmod_mat)
    kroneckerproduct(A::fq_nmod_mat, B::fq_nmod_mat)

Return the Kronecker product of the two matrices `A` and `B`.
"""
⊗(A::fq_nmod_mat, B::fq_nmod_mat) = kronecker_product(A, B)
kron(A::fq_nmod_mat, B::fq_nmod_mat) = kronecker_product(A, B)
tensorproduct(A::fq_nmod_mat, B::fq_nmod_mat) = kronecker_product(A, B)
kroneckerproduct(A::fq_nmod_mat, B::fq_nmod_mat) = kronecker_product(A, B)
# nrows(A::T) where T = size(A, 1)
# ncols(A::T) where T = size(A, 2)

# I think we should avoid length checking here and return it for entire matrix if given
# Hammingweight(v::T) where T <: Union{fq_nmod_mat, gfp_mat, Vector{S}} where S <: Int = count(i->(i != 0), v)
"""
    Hammingweight(v::T) where T <: Union{fq_nmod_mat, Vector{S}} where S <: Int
    weight(v::T) where T <: Union{fq_nmod_mat, Vector{S}} where S <: Int
    wt(v::T) where T <: Union{fq_nmod_mat, Vector{S}} where S <: Int

Return the Hamming weight of `v`.
"""
function Hammingweight(v::T) where T <: Union{fq_nmod_mat, Vector{fq_nmod}, Vector{S}, Adjoint{Int, Vector{Int}}} where S <: Int
    count = 0
    for i in 1:length(v)
        if !iszero(v[i])
            count += 1
        end
    end
    return count
end
weight(v::T) where T <: Union{fq_nmod_mat, Vector{fq_nmod}, Vector{S}, Adjoint{Int, Vector{Int}}} where S <: Int = Hammingweight(v)
wt(v::T) where T <: Union{fq_nmod_mat, Vector{fq_nmod}, Vector{S}, Adjoint{Int, Vector{Int}}} where S <: Int = Hammingweight(v)

function Hammingweight(v::Matrix{Int})
    return sum(v)
end

# TODO: should do the full row, col double loop
function wt(v::Matrix{Int})
    count = 0
    for i in 1:length(v)
        if !iszero(v[1, i])
            count += 1
        end
    end
    return count
end

"""
    wt(f::fq_nmod_poly)

Return the number of nonzero coefficients of the polynomial `f`.
"""
wt(f::fq_nmod_poly) = Hammingweight(collect(coefficients(f)))

"""
    _minwtrow(A::fq_nmod_mat)

Return the minimum weight and corresponding index of the rows of `A`.
"""
function _minwtrow(A::Union{fq_nmod_mat, Matrix{Int}, LinearAlgebra.Adjoint{Int, Matrix{Int}}})
    nr, nc = size(A)
    w = nc
    i = 0
    for r in 1:nr
        wloc = 0
        for c in 1:nc
            iszero(A[r, c]) || (wloc += 1;)
        end
        wloc < w && (w = wloc; i = r;)
    end
    return w, i
end

function _minwtcol(A::LinearAlgebra.Adjoint{Int, Matrix{Int}})
    nr, nc = size(A)
    w = nr
    i = 0
    for c in 1:nc
        wloc = 0
        for r in 1:nr
            iszero(A[r, c]) || (wloc += 1;)
        end
        wloc < w && (w = wloc; i = c;)
    end
    return w, i
end

"""
    Hammingdistance(u::T, v::T) where T <: Union{fq_nmod_mat, Vector{S}} where S <: Int
    distance(u::T, v::T) where T <: Union{fq_nmod_mat, Vector{S}} where S <: Int
    dist(u::T, v::T) where T <: Union{fq_nmod_mat, Vector{S}} where S <: Int

Return the Hamming distance between `u` and `v`.
"""
Hammingdistance(u::T, v::T) where T <: Union{fq_nmod_mat, Vector{S}} where S <: Int = Hammingweight(u .- v)
distance(u::T, v::T) where T <: Union{fq_nmod_mat, Vector{S}} where S <: Int = Hammingweight(u .- v)
dist(u::T, v::T) where T <: Union{fq_nmod_mat, Vector{S}} where S <: Int = Hammingweight(u .- v)

"""
    symplecticinnerproduct(u::fq_nmod_mat, v::fq_nmod_mat)

Return the symplectic inner product of `u` and `v`.
"""
function symplecticinnerproduct(u::fq_nmod_mat, v::fq_nmod_mat)
    (nrows(u) == 1 || ncols(u) == 1) || error("First argument of symplectic inner product is not a vector: dims = $(size(u, 1)).")
    (nrows(v) == 1 || ncols(v) == 1) || error("Second argument of symplectic inner product is not a vector: dims = $(size(v, 1)).")
    length(u) == length(v) || error("Vectors must be the same length in symplectic inner product.")
    iseven(length(u)) || error("Vectors must have even length in symplectic inner product.")
    base_ring(u) == base_ring(v) || error("Vectors must be over the same field in symplectic inner product.")
    ncols = div(length(u), 2)
    return sum([u[i + ncols] * v[i] - v[i + ncols] * u[i] for i in 1:ncols])
end

"""
    aresymplecticorthogonal(A::fq_nmod_mat, B::fq_nmod_mat, symp::Bool=false)

Return `true` if the rows of the matrices `A` and `B` are symplectic orthogonal.

If the optional parameter `symp` is set to `true`, `A` and `B` are assumed to be
in symplectic form over the base field.
"""
function aresymplecticorthogonal(A::fq_nmod_mat, B::fq_nmod_mat, symp::Bool=false)
    E = base_ring(A)
    E == base_ring(B) || error("Matrices in product must both be over the same base ring.")
    if symp
        iseven(ncols(A)) || error("Expected a symplectic input but the first input matrix has an odd number of columns.")
        iseven(ncols(B)) || error("Expected a symplectic input but the second input matrix has an odd number of columns.")
    else
        iseven(degree(E)) || error("The base ring of the given matrices are not a quadratic extension.")
        A = quadratictosymplectic(A)
        B = quadratictosymplectic(B)
    end

    AEuc = hcat(A[:, div(ncols(A), 2) + 1:end], -A[:, 1:div(ncols(A), 2)])
    iszero(AEuc * transpose(B)) || return false
    return true
end

# function traceinnerproduct(u::fq_nmod_mat, v::fq_nmod_mat)
#
# end

"""
    Hermitianinnerproduct(u::fq_nmod_mat, v::fq_nmod_mat)

Return the Hermitian inner product of `u` and `v`.
"""
function Hermitianinnerproduct(u::fq_nmod_mat, v::fq_nmod_mat)
    (nrows(u) == 1 || ncols(u) == 1) || error("First argument of Hermitian inner product is not a vector: dims = $(size(u, 1)).")
    (nrows(v) == 1 || ncols(v) == 1) || error("Second argument of Hermitian inner product is not a vector: dims = $(size(v, 1)).")
    length(u) == length(v) || error("Vectors must be the same length in Hermitian inner product.")
    base_ring(u) == base_ring(v) || error("Vectors must be over the same field in Hermitian inner product.")
    q2 = order(base_ring(u))
    issquare(q2) || error("The Hermitian inner product is only defined over quadratic field extensions.")
    q = Int(sqrt(q2))
    return sum([u[i] * v[i]^q for i in 1:length(u)])
end

"""
    Hermitianconjugatematrix(A::fq_nmod_mat)

Return the Hermitian conjugate of the matrix `A`.
"""
function Hermitianconjugatematrix(A::fq_nmod_mat)
    B = copy(A)
    q2 = order(base_ring(A))
    issquare(q2) || error("The Hermitian conjugate is only defined over quadratic field extensions.")
    q = Int(sqrt(q2))
    return B .^ q
end

"""
    entropy(x::Real)

Return the entropy of the real number `x`.
"""
function entropy(x::Real)
    x != 0 || return 0
    (0 < x <= 1 - 1 / q) || error("Number should be in the range [0, 1 - 1/order(field)].")
    F = parent(x)
    q = order(F)
    return x * (log(q, q - 1) - log(q, x)) - (1 - x) * log(q, 1 - x)
end

"""
    FpmattoJulia(M::fq_nmod_mat)

Return the `fq_nmod_mat` matrix `M` as a Julia Int matrix.
"""
function FpmattoJulia(M::fq_nmod_mat)
    degree(base_ring(M)) == 1 || error("Cannot promote higher order elements to the Ints.")
    A = zeros(Int, size(M))
    for c in 1:ncols(M)
        for r in 1:nrows(M)
            A[r, c] = coeff(M[r, c], 0)
        end
    end
    return A
end

"""
    istriorthogonal(G::fq_nmod_mat, verbose::Bool=false)
    istriorthogonal(G::Matrix{Int}, verbose::Bool=false)

Return `true` if the binary matrix `G` is triorthogonal (modulo 2).

If the optional parameter `verbos` is set to `true`, the first pair or triple of
non-orthogonal rows will be identified on the console.
"""
function istriorthogonal(G::fq_nmod_mat, verbose::Bool=false)
    Int(order(base_ring(G))) == 2 || error("Triothogonality is only defined over 𝔽_2.")
    nr, nc = size(G)
    for r1 in 1:nr
        for r2 in 1:nr
            @views g1 = G[r1, :]
            @views g2 = G[r2, :]
            @views if !iszero(sum([g1[1, i] * g2[1, i] for i in 1:nc]))
                verbose && println("Rows $r1 and $r2 are not orthogonal.")
                return false
            end
        end
    end

    for r1 in 1:nr
        for r2 in 1:nr
            for r3 in 1:nr
                @views g1 = G[r1, :]
                @views g2 = G[r2, :]
                @views g3 = G[r3, :]
                @views if !iszero(sum([g1[1, i] * g2[1, i] * g3[1, i] for i in 1:nc]))
                    verbose && println("Rows $r1, $r2, and $r3 are not orthogonal.")
                    return false
                end
            end
        end
    end
    return true
end

function istriorthogonal(G::Matrix{Int}, verbose::Bool=false)
    nr, nc = size(G)
    for r1 in 1:nr
        for r2 in 1:nr
            @views g1 = G[r1, :]
            @views g2 = G[r2, :]
            @views if !iszero(sum([g1[1, i] * g2[1, i] for i in 1:nc]) % 2)
                verbose && println("Rows $r1 and $r2 are not orthogonal.")
                return false
            end
        end
    end

    for r1 in 1:nr
        for r2 in 1:nr
            for r3 in 1:nr
                @views g1 = G[r1, :]
                @views g2 = G[r2, :]
                @views g3 = G[r3, :]
                @views if !iszero(sum([g1[1, i] * g2[1, i] * g3[1, i] for i in 1:nc]) % 2)
                    verbose && println("Rows $r1, $r2, and $r3 are not orthogonal.")
                    return false
                end
            end
        end
    end
    return true
end

function _quotientspace(big::fq_nmod_mat, small::fq_nmod_mat)
    F = base_ring(big)
    V = VectorSpace(F, ncols(big))
    U, UtoV = sub(V, [V(small[i, :]) for i in 1:nrows(small)])
    W, WtoV = sub(V, [V(big[i, :]) for i in 1:nrows(big)])
    gensofUinW = [preimage(WtoV, UtoV(g)) for g in gens(U)]
    UinW, _ = sub(W, gensofUinW)
    Q, WtoQ = quo(W, UinW)
    iszero(dim(Q)) && (return zero_matrix(F, 1, ncols(big));)
    C2modC1basis = [WtoV(x) for x in [preimage(WtoQ, g) for g in gens(Q)]]
    Fbasis = [[F(C2modC1basis[j][i]) for i in 1:AbstractAlgebra.dim(parent(C2modC1basis[1]))] for j in 1:length(C2modC1basis)]
    return matrix(F, length(Fbasis), length(Fbasis[1]), vcat(Fbasis...))
end

#############################
  # Quantum Helper Functions
#############################

function printstringarray(A::Vector{String}, withoutIs=false)
    for a in A
        if !withoutIs
            println(a)
        else
            for i in a
                if i == 'I'
                    print(' ')
                else
                    print(i)
                end
            end
            print('\n')
        end
    end
end
printchararray(A::Vector{Vector{Char}}, withoutIs=false) = printstringarray(setchartostringarray(A), withoutIs)
printsymplecticarray(A::Vector{Vector{T}}, withoutIs=false) where T <: Int = printstringarray(setsymplectictostringarray(A), withoutIs)

"""
    pseudoinverse(M::fq_nmod_mat)

Return the pseudoinverse of a stabilizer matrix `M` over a quadratic extension.

Note that this is not the Penrose-Moore pseudoinverse.
"""
function pseudoinverse(M::fq_nmod_mat)
    # let this fail elsewhere if not actually over a quadratic extension
    if degree(base_ring(M)) != 1
        M = transpose(quadratictosymplectic(M))
    else
        M = transpose(M)
    end

    nr, nc = size(M)
    MS = MatrixSpace(base_ring(M), nr, nr)
    _, E = rref(hcat(M, MS(1)))
    E = E[:, (nc + 1):end]
    pinv = E[1:nc, :]
    dual = E[nc + 1:nr, :]

    # verify
    _, Mrref = rref(M)
    MScols = MatrixSpace(base_ring(M), nc, nc)
    E * M == Mrref || error("Pseudoinverse calculation failed (transformation incorrect).")
    Mrref[1:nc, 1:nc] == MScols(1) || error("Pseudoinverse calculation failed (failed to get I).")
    iszero(Mrref[nc + 1:nr, :]) || error("Pseudoinverse calculation failed (failed to get zero).")
    pinv * M == MScols(1) || error("Pseudoinverse calculation failed (eq 1).")
    transpose(M) * transpose(pinv) == MScols(1) || error("Pseudoinverse calculation failed (eq 2).")
    iszero(transpose(M) * transpose(dual)) || error("Failed to correctly compute dual (rhs).")
    iszero(dual * M) || error("Failed to correctly compute dual (lhs).")
    return pinv
end

"""
    quadratictosymplectic(M::fq_nmod_mat)

Return the matrix `M` converted from the quadratic to the symplectic form.
"""
function quadratictosymplectic(M::fq_nmod_mat)
    E = base_ring(M)
    iseven(degree(E)) || error("The base ring of the given matrix is not a quadratic extension.")
    F, _ = FiniteField(Int(characteristic(E)), div(degree(E), 2), "ω")
    nr = nrows(M)
    nc = ncols(M)
    Msym = zero_matrix(F, nr, 2 * nc)
    for c in 1:nc
        for r in 1:nr
            # TODO: benchmark this without the branching
            if !iszero(M[r, c])
                Msym[r, c] = F(coeff(M[r, c], 0))
                Msym[r, c + nc] = F(coeff(M[r, c], 1))
            end
        end
    end
    return Msym
end

"""
    symplectictoquadratic(M::fq_nmod_mat)

Return the matrix `M` converted from the symplectic to the quadratic form.
"""
function symplectictoquadratic(M::fq_nmod_mat)
    iseven(ncols(M)) || error("Input to symplectictoquadratic is not of even length.")
    nr = nrows(M)
    nc = div(ncols(M), 2)
    F = base_ring(M)
    E, ω = FiniteField(Int(characteristic(F)), 2 * degree(F), "ω")
    ϕ = embed(F, E)
    Mquad = zero_matrix(E, nr, nc)
    for c in 1:nc
        for r in 1:nr
            Mquad[r, c] = ϕ(M[r, c]) + ϕ(M[r, c + nc]) * ω
        end
    end
    return Mquad
end

function _Paulistringtosymplectic(str::T) where T <: Union{String, Vector{Char}}
    n = length(str)
    F, _ = FiniteField(2, 1, "ω")
    sym = zero_matrix(F, 1, 2 * n)
    for (i, c) in enumerate(str)
        if c == 'X'
            sym[1, i] = F(1)
        elseif c == 'Z'
            sym[1, i + n] = F(1)
        elseif c == 'Y'
            sym[1, i] = 1
            sym[1, i + n] = F(1)
        elseif c != 'I'
            error("Encountered non-{I, X, Y, Z} character in Pauli string. This function is only defined for binary strings.")
        end
    end
    return sym
end
_Paulistringtosymplectic(A::Vector{T}) where T <: Union{String, Vector{Char}} = vcat([_Paulistringtosymplectic(s) for s in A]...)
_Paulistringstofield(str::T) where T <: Union{String, Vector{Char}} = symplectictoquadratic(_Paulistringtosymplectic(str))
_Paulistringstofield(A::Vector{T}) where T <: Union{String, Vector{Char}} = vcat([_Paulistringstofield(s) for s in A]...)
# need symplectictoPaulistring
# quadratictoPaulistring

# charvec::Union{Vector{nmod}, Missing}=missing)
function _processstrings(SPauli::Vector{T}) where T <: Union{String, Vector{Char}}
    # Paulisigns = Vector{Int}()
    StrPaulistripped = Vector{String}()
    for (i, s) in enumerate(SPauli)
        if s[1] ∈ ['I', 'X', 'Y', 'Z']
            # append!(Paulisigns, 1)
            push!(StrPaulistripped, s)
        elseif s[1] == '+'
            # append!(Paulisigns, 1)
            push!(StrPaulistripped, s[2:end])
        elseif s[1] == '-'
            # append!(Paulisigns, -1)
            push!(StrPaulistripped, s[2:end])
        else
            error("The first element of Pauli string $i is neither a Pauli character or +/-: $s.")
        end
    end

    n = length(StrPaulistripped[1])
    for s in StrPaulistripped
        for i in s
            i ∈ ['I', 'X', 'Y', 'Z'] || error("Element of provided Pauli string is not a Pauli character: $s.")
        end
        length(s) == n || error("Not all Pauli strings are the same length.")
    end

    # if !ismissing(charvec)
    #     2 * n == length(charvec) || error("The characteristic value is of incorrect length.")
    #     R = ResidueRing(Nemo.ZZ, 4)
    #     for s in charvec
    #         modulus(s) == modulus(R) || error("Phases are not in the correct ring.")
    #     end
    # else
    #     R = ResidueRing(Nemo.ZZ, 4)
    #     charvec = [R(0) for _ in 1:2 * n]
    # end
    return StrPaulistripped#, charvec
end

function largestconsecrun(arr::Vector{Int})
    n = length(arr)
    maxlen = 1
    for i = 1:n
        mn = arr[i]
        mx = arr[i]

        for j = (i + 1):n
            mn = min(mn, arr[j])
            mx = max(mx, arr[j])

            if (mx - mn) == (j - i)
                maxlen = max(maxlen, mx - mn + 1)
            end
        end
    end

    return maxlen
end

function _removeempty(A::fq_nmod_mat, type::String)
    type ∈ ["rows", "cols"] || error("Unknown type in _removeempty; expected: `rows` or `cols`, received: $type")
    del = Vector{Int}()
    if type == "rows"
        for r in 1:nrows(A)
            if iszero(A[r, :])
                append!(del, r)
            end
        end
        return isempty(del) ? A : A[setdiff(1:nrows(A), del), :]
    else
        for c in 1:ncols(A)
            if iszero(A[:, c])
                append!(del, c)
            end
        end
        return isempty(del) ? A : A[:, setdiff(1:ncols(A), del)]
    end
end

function _rref_no_col_swap(M::fq_nmod_mat, rowrange::UnitRange{Int}, colrange::UnitRange{Int})
    isempty(rowrange) && throw(ArgumentError("The row range cannot be empty in _rref_no_col_swap."))
    isempty(colrange) && throw(ArgumentError("The column range cannot be empty in _rref_no_col_swap."))
    A = deepcopy(M)

    i = rowrange.start
    j = colrange.start
    nr = rowrange.stop
    nc = colrange.stop
    while i <= nr && j <= nc
        # find first pivot
        ind = 0
        for k in i:nr
            if !iszero(A[k, j])
                ind = k
                break
            end
        end

        if !iszero(ind)
            # normalize pivot
            if !isone(A[ind, j])
                A[ind, :] *= inv(A[ind, j])
            end

            # swap to put the pivot in the next row
            if ind != i
                A[i, :], A[ind, :] = A[ind, :], A[i, :]
            end

            # eliminate
            for k = rowrange.start:nr
                if k != i
                    # do a manual loop here to reduce allocations
                    d = A[k, j]
                    @simd for l = 1:nc
                        A[k, l] = (A[k, l] - d * A[i, l])
                    end
                end
            end
            i += 1
        end
        j += 1
    end
    return A
end

function _rref_col_swap(M::fq_nmod_mat, rowrange::UnitRange{Int}, colrange::UnitRange{Int})
    isempty(rowrange) && throw(ArgumentError("The row range cannot be empty in _rref_no_col_swap."))
    isempty(colrange) && throw(ArgumentError("The column range cannot be empty in _rref_no_col_swap."))
    A = deepcopy(M)
    # permutation matrix required to return to rowspace if column swap done
    P = missing
    ncA = ncols(A)

    rnk = 0
    i = rowrange.start
    j = colrange.start
    nr = rowrange.stop
    nc = colrange.stop
    while i <= nr && j <= nc
        # find first pivot
        ind = 0
        for k in i:nr
            if !iszero(A[k, j])
                ind = k
                break
            end
        end

        # need to column swap
        if iszero(ind)
            for k in j + 1:nc
                for l in i:nr
                    if !iszero(A[l, k])
                        ismissing(P) && (P = identity_matrix(base_ring(A), nc);)
                        swap_cols!(A, k, j)
                        swap_rows!(P, k, j)
                        ind = l
                        break
                    end
                end
            end
        end

        # if true, the rest of the submatrix is zero
        if iszero(ind)
            return rnk, A, P
        else
            # normalize pivot
            if !isone(A[ind, j])
                A[ind, :] *= inv(A[ind, j])
            end

            # swap to put the pivot in the next row
            ind != i && swap_rows!(A, ind, i)

            # eliminate
            for k = rowrange.start:nr
                if k != i
                    # do a manual loop here to reduce allocations
                    d = A[k, j]
                    @simd for l = 1:ncA
                        A[k, l] = (A[k, l] - d * A[i, l])
                    end
                end
            end
        end
        i += 1
        j += 1
        rnk += 1
    end
    return rnk, A, P
end

function digitstoint(x::Vector{Int}, base::Int=2)
    res = 0
    lenx = length(x)
    for i in 1:lenx
        res = x[i] + base * res
    end
    return res
end

"""
    polytocircmatrix(f::AbstractAlgebra.Generic.Res{fq_nmod_poly})

Return the circulant matrix whose first column is the coefficients of `f`.
"""
function polytocircmatrix(f::AbstractAlgebra.Generic.Res{fq_nmod_poly})
    R = parent(f)
    S = base_ring(R)
    F = base_ring(S)
    g = modulus(R)
    l = degree(g)
    g == gen(S)^l - 1 || throw(ArgumentError("Residue ring not of the form x^l - 1."))
    # gcd(l, Int(characteristic(F))) == 1 || throw(ArgumentError("Residue ring over F_q[x] must have modulus x^l - 1 with gcd(l, q) = 1."))

    A = zero_matrix(F, l, l)
    fcoeffs = zero_matrix(F, l, 1)
    temp = collect(coefficients(Nemo.lift(f)))
    fcoeffs[1:length(temp), 1] = temp
    A[:, 1] = fcoeffs
    for c in 2:l
        A[:, c] = circshift(fcoeffs, c - 1)
    end
    return A
end

"""
    lift(A::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}})

Return the matrix whose polynomial elements are converted to circulant matrices over the base field.
"""
function lift(A::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}})
    R = parent(A[1, 1])
    S = base_ring(R)
    F = base_ring(S)
    g = modulus(R)
    l = degree(g)
    g == gen(S)^l - 1 || throw(ArgumentError("Residue ring not of the form x^l - 1."))
    # gcd(l, Int(characteristic(F))) == 1 || throw(ArgumentError("Residue ring over F_q[x] must have modulus x^l - 1 with gcd(l, q) = 1."))

    nr, nc = size(A)
    Alift = zero_matrix(F, nr * l, nc * l)
    for c in 1:nc
        for r in 1:nr
            if !iszero(A[r, c])
                Alift[(r - 1) * l + 1:r * l, (c - 1) * l + 1:c * l] = polytocircmatrix(A[r, c])
            end
        end
    end
    return Alift
end

#############################
        # Finite Fields
#############################

"""
    tr(x::fq_nmod, K::FqNmodFiniteField, verify::Bool=false)

Return the relative trace of `x` from its base field to the field `K`.

If the optional parameter `verify` is set to `true`, the two fields are checked
for compatibility.
"""
function tr(x::fq_nmod, K::FqNmodFiniteField, verify::Bool=false)
    L = parent(x)
    q = order(K)
    if verify
        # # shouldn't need Int casting here but just in case...
        # Int(characteristic(L)) == Int(characteristic(K)) || error("The given field is not a subfield of the base ring of the element.")
        # degree(L) % degree(K) == 0 || error("The given field is not a subfield of the base ring of the element.")
        flag, m = isextension(L, K)
        flag || throw(ArgumentError("The given field is not a subfield of the base ring of the matrix."))
    end
    n = div(degree(L), degree(K))
    return sum([x^(q^i) for i in 0:(n - 1)])
end

function _expandelement(x::fq_nmod, K::FqNmodFiniteField, basis::Vector{fq_nmod}, verify::Bool=false)
    return [tr(x * i) for i in basis] #, K, verify
end

function _expandrow(row::fq_nmod_mat, K::FqNmodFiniteField, basis::Vector{fq_nmod}, verify::Bool=false)
    new_row = _expandelement(row[1], K, basis, verify)
    for i in 2:ncols(row)
        new_row = vcat(new_row, _expandelement(row[i], K, basis, verify))
    end
    return matrix(K, 1, length(new_row), new_row)
end

"""
    expandmatrix(M::fq_nmod_mat, K::FqNmodFiniteField, basis::Vector{fq_nmod})

Return the matrix constructed by expanding the elements of `M` to the subfield
`K` using the provided `basis` for the base ring of `M` over `K`.
"""
function expandmatrix(M::fq_nmod_mat, K::FqNmodFiniteField, basis::Vector{fq_nmod})
    L = base_ring(M)
    L == K && return M
    Int(characteristic(L)) == Int(characteristic(K)) || error("The given field is not a subfield of the base ring of the element.")
    degree(L) % degree(K) == 0 || error("The given field is not a subfield of the base ring of the element.")
    n = div(degree(L), degree(K))
    n == length(basis) || error("Provided basis is of incorrect size for the given field and subfield.")
    # should really check if it is a basis
    flag, m = isextension(L, K)
    flag || throw(ArgumentError("The given field is not a subfield of the base ring of the matrix."))
    m == length(basis) || throw(ArgumentError("Basis does not have length degree of the extension."))
    flag, _ = _isbasis(L, basis, Int(order(K)))
    flag || throw(ArgumentError("The provided vector is not a basis for the extension."))
    return vcat([_expandrow(M[r, :], K, basis) for r in 1:nrows(M)]...)
end

"""
    quadraticresidues(q::Int, n::Int)

Return the set of quadratic resides and quadratic non-residues of `q` and `n`.
"""
# TODO: improve this description
function quadraticresidues(q::Int, n::Int)
    isodd(n) && isprime(n) || error("n must be an odd prime in quadratic residues")
    q^div(n - 1, 2) % n == 1 || error("q^(n - 1)/2 ≅ 1 mod n in quadratic residues")

    # F, _ = FiniteField(n, 1, "α")
    # elms = collect(F)
    # # skip 0
    # qres = unique!([i^2 for i in elms[2:end]])
    # nqres = setdiff!(elms[2:end], qres)
    # return qres, nqres

    # don't want this returning in the field
    qres = sort!(unique!([i^2 % n for i in 1:n - 1]))
    nqres = setdiff(1:(n - 1), qres)
    return qres, nqres
end

function _isbasis(E::FqNmodFiniteField, basis::Vector{fq_nmod}, q::Int)
    m = length(basis)
    M = MatrixSpace(E, m, m)
    B = M(0)
    for r in 1:m
        for c in 1:m
            B[r, c] = basis[r]^(q^(c - 1))
        end
    end
    iszero(det(B)) && return false, missing
    
    try
        Binv = inv(B)
        λ = [Binv[1, i] for i in 1:m]
        return true, λ
    catch
        return false, missing
    end
end

"""
    isextension(E::FqNmodFiniteField, F::FqNmodFiniteField)

Return `true` if `E/F` is a valid field extension.
"""
function isextension(E::FqNmodFiniteField, F::FqNmodFiniteField)
    p = Int(characteristic(E))
    Int(characteristic(F)) == p || return false, missing
    degE = degree(E)
    degF = degree(F)
    degE % degF == 0 || return false, missing
    return true, div(degE, degF)
    # the below allows you to embed GF(2) into GF(5) without error but is not an extension
    # try
    #     embed(F, E)
    #     return true, div(degree(E), degree(F))
    # catch
    #     return false, missing
    # end
end

"""
    isbasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod})

Return `true` and the dual (complementary) basis if `basis` is a basis for `E/F`,
otherwise return `false, missing`.
"""
function isbasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod})
    flag, m = isextension(E, F)
    flag || throw(ArgumentError("Second field is not a subfield of the first."))
    length(basis) == m || throw(ArgumentError("Basis does not have length degree of the extension."))
    for i in 1:m
        parent(basis[i]) == E || throw(ArgumentError("The basis must be elements of the extension field."))
    end
    return _isbasis(E, basis, Int(order(F)))
end

"""
    primitivebasis(E::FqNmodFiniteField, F::FqNmodFiniteField)

Return a primitive basis for `E/F` and its dual (complementary) basis.
"""
function primitivebasis(E::FqNmodFiniteField, F::FqNmodFiniteField)
    flag, m = isextension(E, F)
    flag || throw(ArgumentError("Second field is not a subfield of the first."))
    α = gen(E)
    basis = [α^i for i in 0:m - 1]
    flag, λ = _isbasis(E, basis, Int(order(F)))
    return basis, λ
end
# these are slightly different
# polynomialbasis(E::FqNmodFiniteField, F::FqNmodFiniteField) = primitivebasis(E, F)
# monomialbasis(E::FqNmodFiniteField, F::FqNmodFiniteField) = primitivebasis(E, F)

"""
    normalbasis(E::FqNmodFiniteField, F::FqNmodFiniteField)

Return a normal basis for `E/F` and its dual (complementary) basis.
"""
# "Normal Bases over Finite Fields" by Shuhong Gao has algorithms for this but they are
# complicated for the field sizes intended in this work
function normalbasis(E::FqNmodFiniteField, F::FqNmodFiniteField)
    flag, m = isextension(E, F)
    flag || throw(ArgumentError("Second field is not a subfield of the first."))

    q = Int(order(F))
    elms = collect(E)
    for e in elms
        basis = [e^(q^i) for i in 0:m - 1]
        flag, dualbasis = _isbasis(E, basis, q)
        flag && return basis, dualbasis
    end
    error("Somehow failed to final a normal element for the extension.")
end

"""
    dualbasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod})
    complementarybasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod})

Return the dual (complentary) basis of `basis` for the extension `E/F`.
"""
function dualbasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod})
    flag, λ = isbasis(E, F, basis)
    flag || throw(ArgumentError("The provided vector is not a basis for the extension."))
    return λ
end
complementarybasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod}) = dualbasis(E, F, basis)

"""
    verifydualbasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod}, dualbasis::Vector{fq_nmod})
    verifycomplementarybasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod}, dualbasis::Vector{fq_nmod})

Return `true` if `basis` is the dual of `dualbasis` for `E/F`, otherwise return `false`.
"""
function verifydualbasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod}, dualbasis::Vector{fq_nmod})
    flag, m = isextension(E, F)
    flag || throw(ArgumentError("Second field is not a subfield of the first."))

    m = length(basis)
    length(dualbasis) == m || throw(ArgumentError("The basis and dual basis must have the same length."))
    E = parent(basis[1])
    for i in 1:m
        parent(basis[i]) == E || throw(ArgumentError("Elements must be over the same field."))
        parent(dualbasis[i]) == E || throw(ArgumentError("Elements must be over the same field."))
    end

    q = Int(order(F))
    M = MatrixSpace(E, m, m)
    B = M(0)
    for r in 1:m
        for c in 1:m
            B[r, c] = basis[r]^(q^(c - 1))
        end
    end
    Binv = M(0)
    for r in 1:m
        for c in 1:m
            Binv[r, c] = dualbasis[c]^(q^(r - 1))
        end
    end
    return B * Binv == M(1)
end
verifycomplementarybasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod}, dualbasis::Vector{fq_nmod}) = verifydualbasis(E, F, basis, dualbasis)

"""
    isequivalentbasis(basis::Vector{fq_nmod}, basis2::Vector{fq_nmod})

Return `true` if `basis` is a scalar multiple of `basis2`.`
"""
function isequivalentbasis(basis::Vector{fq_nmod}, basis2::Vector{fq_nmod})
    m = length(basis)
    length(basis2) == m || throw(ArgumentError("The two vectors must have the same length."))
    c = basis[1] * basis2[1]^-1
    E = parent(basis[1])
    for i in 1:m
        parent(basis[i]) == E || throw(ArgumentError("Elements must be over the same field."))
        parent(basis2[i]) == E || throw(ArgumentError("Elements must be over the same field."))
        # for logical consistency should probably do this in a separate loop
        basis[i] == c * basis2[i] || return false
    end
    return true
end

"""
    isselfdualbasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod})

Return `true` if `basis` is equal to its dual.
"""
function isselfdualbasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod})
    flag, λ = isbasis(E, F, basis)
    flag || throw(ArgumentError("The provided vector is not a basis for the extension."))
    return basis == λ
end

"""
    isprimitivebasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod})

Return `true` if `basis` is a primitive basis for `E/F`.
"""
function isprimitivebasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod})
    flag, _ = isbasis(E, F, basis)
    flag || throw(ArgumentError("The provided vector is not a basis for the extension."))
    isone(basis[1]) ? (x = basis[2];) : (x = basis[1];)
    m = length(basis)
    for i in 0:m - 1
        x^i ∈ basis || return false
    end
    return true
end

"""
    isnormalbasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod})

Return `true` if `basis` is a normal basis for `E/F`.
"""
function isnormalbasis(E::FqNmodFiniteField, F::FqNmodFiniteField, basis::Vector{fq_nmod})
    flag, _ = isbasis(E, F, basis)
    flag || throw(ArgumentError("The provided vector is not a basis for the extension."))
    isone(basis[1]) ? (return false;) : (x = basis[1];)
    m = length(basis)
    q = Int(order(F))
    for i in 0:m - 1
        x^(q^i) ∈ basis || return false
    end
    return true
end

# The finite field F_q^n has a pair of self-dual bases for the following parameters.
# (1) q is an even prime power.
# (2) q is an odd prime power and n = 2k + 1.
# (Imamura 1983) The finite field F_q^n has no self-dual power bases.

#############################
          # Graphs
#############################

"""
    isregular(G::SimpleGraph{Int})

Return `true` if `G` is regular.
"""
function isregular(G::SimpleGraph{Int})
    fadjlist = G.fadjlist
    deg = length(fadjlist[1])
    for v in fadjlist
        length(v) == deg || return false
    end
    return true
end

"""
    edgevertexincidencematrix(G::SimpleGraph{Int})

Return the edge-vertex incidence matrix of `G` along with the vertex incides of the left
and right bipartition.
"""
function edgevertexincidencematrix(G::SimpleGraph{Int})
    I = Array(incidence_matrix(G))
    nr, nc = size(I)
    Itr = transpose(I)
    B = vcat(hcat(zeros(Int, nc, nc), Itr), hcat(I, zeros(Int, nr, nr)))
    return B, collect(1:nc), collect(nc + 1:nr + nc)
end

"""
    edgevertexincidencegraph(G::SimpleGraph{Int})

Return the edge-vertex incidence graph of `G` along with the vertex incides of the left
and right bipartition.
"""
function edgevertexincidencegraph(G::SimpleGraph{Int})
    B, left, right = edgevertexincidencematrix(G)
    return SimpleGraph(B), left, right
end

"""
    isvalidbipartition(G::SimpleGraph{Int}, left::Vector{Int}, right::Vector{Int})

Return `true` if the vertices indexed by `left` and `right` form a valid bipartition for `G`.
"""
function isvalidbipartition(G::SimpleGraph{Int}, left::Vector{Int}, right::Vector{Int})
    nv(G) == length(left) + length(right) || throw(ArgumentError("Too few vertices in lists."))
    l = sort(left)
    r = sort(right)
    temp = l ∩ r # can manually do using sorted knowledge if this is slow
    isempty(temp) || throw(ArgumentError("Bipartition must be disjoint."))
    # could consider getting rid of the above in exchange for changing the elseif to an if

    for (i, v) in enumerate(G.fadjlist)
        if insorted(i, l)
            for e in v
                insorted(e, r) || return false
            end
        elseif insorted(i, r)
            for e in v
                insorted(e, l) || return false
            end
        else
            return false
        end
    end
    return true
end

"""
    extractbipartition(G::SimpleGraph{Int})

Return two vectors representing the vertex indices of each side of the bipartition.
"""
function extractbipartition(G::SimpleGraph{Int})
    temp = bipartite_map(G)
    # this is the definition of the function is_bipartite in Graphs.jl
    length(temp) == nv(G) || throw(ArgumentError("Input graph is not bipartite."))
    left = Vector{Int}()
    right = Vector{Int}()
    for i in temp
        i == 1 ? (push!(i, left);) : (push!(i, right);)
    end
    return left, right
end









# #=
# Example of using the repeated iterator inside of product.
#
# It turns out that this is faster than the Nemo iterator and doesn't allocate.
#
# julia> @benchmark for i in Base.Iterators.product(Base.Iterators.repeated(0:1, 10)...) i end
# BenchmarkTools.Trial: 10000 samples with 137 evaluations.
#  Range (min … max):  713.022 ns …  1.064 μs  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     755.949 ns              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   760.380 ns ± 24.121 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%
#
#  Memory estimate: 0 bytes, allocs estimate: 0.
#
# julia> @benchmark for i in Nemo.AbstractAlgebra.ProductIterator([0:1 for _ in 1:10]) i end
# BenchmarkTools.Trial: 10000 samples with 1 evaluation.
#  Range (min … max):  34.064 μs …   2.604 ms  ┊ GC (min … max):  0.00% … 97.51%
#  Time  (median):     36.970 μs               ┊ GC (median):     0.00%
#  Time  (mean ± σ):   46.342 μs ± 124.916 μs  ┊ GC (mean ± σ):  16.57% ±  6.04%
#
#  Memory estimate: 176.50 KiB, allocs estimate: 2051.
#
# julia> @benchmark for i in Base.Iterators.product([0:1 for _ in 1:10]...) i end
# BenchmarkTools.Trial: 10000 samples with 1 evaluation.
#  Range (min … max):  53.741 μs …   1.465 ms  ┊ GC (min … max):  0.00% … 87.86%
#  Time  (median):     63.790 μs               ┊ GC (median):     0.00%
#  Time  (mean ± σ):   76.919 μs ± 104.655 μs  ┊ GC (mean ± σ):  12.40% ±  8.59%
#
#  Memory estimate: 432.88 KiB, allocs estimate: 2061.
# =#
#
#
# # Gray code iterator, naive formula, gives Ints instead of vectors
#
# struct GrayCodeNaive
#     n::Int
# end
#
# Base.iterate(G::GrayCodeNaive) = G.n < 64 ? (0, 1) : error("Don't handle cases this large")
#
# function Base.iterate(G::GrayCodeNaive, k)
#     k == 2^G.n && return nothing
#     return (k ⊻ (k >> 1), k + 1)
# end
#
# Base.length(G::GrayCodeNaive) = 2^G.n
#
# #=
# Benchmark result:
#
# julia> @benchmark for g in GrayCodeNaive(25) g end
# BenchmarkTools.Trial: 25 samples with 1 evaluation.
#  Range (min … max):  200.014 ms … 202.460 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     200.305 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   200.626 ms ± 659.369 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
#
#  Memory estimate: 0 bytes, allocs estimate: 0.
# =#
#
#
#
# # Gray code iterator, chooses next based on previous, gives Ints instead of vectors
#

# #=
# Benchmark result:
#
# julia> @benchmark for g in GrayCode(25) g end
# BenchmarkTools.Trial: 15 samples with 1 evaluation.
#  Range (min … max):  348.661 ms … 349.411 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     349.050 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   349.048 ms ± 225.710 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
#
#  Memory estimate: 0 bytes, allocs estimate: 0.
# =#
