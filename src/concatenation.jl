# Copyright (c) 2023 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # Classical
#############################

#############################
        # constructors
#############################

# TODO: give control over expansion basis
"""
    ∘(Cout::AbstractLinearCode, Cin::AbstractLinearCode)
    concatenate(Cout::AbstractLinearCode, Cin::AbstractLinearCode)

Return the concatenation of `Cout` and `Cin`.
"""
function concatenate(Cout::AbstractLinearCode, Cin::AbstractLinearCode)
    Fout = Cout.F
    Fin = Cin.F
    β, λ = missing, missing
    if Int(order(Fout)) != Int(order(Fin))
        flag, deg = isextension(Fout, Fin)
        flag || throw(ArgumentError("Galois concatenation requires the outer code to be over an extension field of the inner code"))
        deg % Cin.k == 0 || Cout.n % Cin.k == 0 || 
            throw(ArgumentError("Inner dimension must divide outer length or extension degree"))
        Gout = generatormatrix(Cout, true)
        ismissing(Cout.Pstand) || (Gout = Gout * Cout.Pstand)
        
        β, λ = primitivebasis(Fout, Fin)
        D = _expansiondict(Fout, Fin, λ)
        Gout = _expandmatrix(Gout, D, deg)
        type = :expanded
    else
        Cout.n % Cin.k == 0 || throw(ArgumentError("Inner dimension must divide outer length"))

        Fout == Fin || (Cout = changefield(Cout, Fin);)
        Gout = generatormatrix(Cout, true)
        ismissing(Cout.Pstand) || (Gout = Gout * Cout.Pstand)
        type = :same
    end
    
    Gin = generatormatrix(Cin, true)
    ismissing(Cin.Pstand) || (Gin = Gin * Cin.Pstand)
    G = _concatenatedgeneratormatrix(Gout, Gin)
    Gstand, Hstand, P, k = _standardform(G)
    H = ismissing(P) ? Hstand : Hstand * P
    ub1, _ = _minwtrow(G)
    ub2, _ = _minwtrow(Gstand)
    ub = min(ub1, ub2)

    C = ConcatenatedCode(Cout, Cin, type, β, λ, Fin, ncols(G), k, missing, 1, ub, G, H, Gstand, Hstand, P, missing)
    # TODO: distance check here
    # for a lower bound on the distance, count the number of pieces of Gout in _concatenatedgeneratormatrix
    # with full rank and multiply by the distance of the inner code

    return C
end
∘(Cout::AbstractLinearCode, Cin::AbstractLinearCode) = concatenate(Cout, Cin)

function concatenate(outers::Vector{T}, inners::Vector{T}) where T <: AbstractLinearCode
    isempty(outers) && throw(ArgumentError("List of codes cannot be empty"))
    length(outers) == length(inners) || throw(ArgumentError("Must have the same number of inner and outer codes"))
    for i in 2:length(inners)
        # inners[i] ⊆ inners[i - 1] || throw(ArgumentError("The inner subcodes must be in a decreasing nested sequence"))
        inners[i - 1] ⊆ inners[i] || throw(ArgumentError("The inner subcodes must be in a decreasing nested sequence"))
    end
    F = first(inners).F
    nin = first(inners).n
    # iszero(generatormatrix(inners[end])) || push!(inners, ZeroCode(F, nin))

    ordF = Int(order(F))
    for (i, Cout) in enumerate(outers)
        if Int(order(Cout.F)) == ordF
            Cout.F != F || (outers[i] = changefield(Cout, F);)
        elseif issubfield(F, Cout.F)
            # TODO: expansion step here
            throw(ArgumentError("Only implemented for pre-expanded outer codes currently"))
        else
            throw(ArgumentError("Cannot connect outer code $i field to inner code field"))
        end
    end

    # Are the outer matrices the right size?
    # Need this check since we can take pre-expanded matrices
    nout = divexact(outers[1].n, inners[1].k)
    for i in 2:length(outers)
        nout == divexact(outers[i].n, inners[i].k - inners[i - 1].k) || throw(ArgumentError("The outer matrices are not of the correct size"))
    end

    B = [generatormatrix(inners[1])]
    for i in 2:length(inners)
        Gi = generatormatrix(inners[i])
        Gim1 = generatormatrix(inners[i - 1])
        push!(B, _quotientspace(Gi, Gim1, :VS))
    end

    G1 = reduce(directsum, generatormatrix(C) for C in outers)
    G2 = zero_matrix(F, ncols(G1), nin * nout)
    z = 1
    for i in eachindex(inners)
        for j in 0:nout - 1
            rows = range(z, z + size(B[i], 1) - 1)
            cols = range(j * nin + 1, (j + 1) * nin)
            G2[rows, cols] = B[i]
            z += nrows(B[i])
        end
    end

    return LinearCode(G1 * G2)
end
generalizedconcatenation(outers::Vector{T}, inners::Vector{T}) where T <: AbstractLinearCode = concatenate(outers, inners)
# BZ and Z, multilevel, cascade?

#############################
      # getter functions
#############################

"""
    innercode(C::AbstractConcatenatedCode)

Return the inner code of the concatenation.
"""
innercode(C::AbstractConcatenatedCode) = C.Cin

"""
    outercode(C::AbstractConcatenatedCode)

Return the outer code of the concatenation.
"""
outercode(C::AbstractConcatenatedCode) = C.Cout

"""
    expansionbasis(C::AbstractConcatenatedCode)

Return the basis used to expanded the outer code, if it exists; otherwise return `missing`.
"""
expansionbasis(C::AbstractConcatenatedCode) = C.basis

"""
    expansiondualbasis(C::AbstractConcatenatedCode)

Return the dual basis used to expanded the outer code, if it exists; otherwise return `missing`.
"""
expansiondualbasis(C::AbstractConcatenatedCode) = C.dualbasis

"""
    concatenationtype(C::AbstractConcatenatedCode)

Return `:expanded`, `:same`, or `:generalized` depending on the type of concatenation.
"""
concatenationtype(C::AbstractConcatenatedCode) = C.type

#############################
      # setter functions
#############################

#############################
     # general functions
#############################

function _concatenatedgeneratormatrix(A::CTMatrixTypes, B::CTMatrixTypes) # where T <: CTMatrixTypes
    nrA, ncA = size(A)
    nrB, ncB = size(B)
    t = div(ncA, nrB)
    M = zero_matrix(base_ring(A), nrA, t * ncB)
    for i in 1:t
        M[:, ncB * (i - 1) + 1: ncB * i] = view(A, :, nrB * (i - 1) + 1:nrB * i) * B
    end
    return M
end

# TODO: untested, little chance this works without error
"""
    encode(C::AbstractConcatenatedCode, v::Union{CTMatrixTypes, Vector{Int}})

Return the encoding of `v` into `C`, where `v` is either a valid input for the outer code or the full code.
"""
function encode(C::AbstractConcatenatedCode, v::Union{CTMatrixTypes, Vector{Int}})
    if typeof(v) <: CTMatrixTypes
        nrv, ncv = size(v)
        nrv == 1 || ncv == 1 || throw(ArgumentError("Vector has incorrect dimension"))
        nrv != 1 && ncv == 1 ? (w = transpose(v);) : (w = v;)
        ncw = ncols(w)
        if ncw == C.Cout.k
            # TODO: should check order and then convert if they are the same but different pointers
            base_ring(w) == C.Cout.F || throw(ArgumentError("Vector must have the same base ring as the outer code."))
            Gout = generatormatrix(C.Cout, true)
            ismissing(C.Cout.Pstand) || (Gout = Gout * C.Cout.Pstand)
            temp = w * Gout
            if C.type == :expanded
                D = _expansiondict(C.Cout.F, C.Cin.F, C.dualbasis)
                temp = _expandmatrix(temp, D, div(degree(C.Cout.F), degree(C.Cin.F)))
                # TODO: this is no longer a vector, only want 1 * temp row instead of full basis?
                temp = temp[1, :]
            end
            # should automatically now be in field of inner code
            Gin = generatormatrix(C.Cin, true)
            ismissing(C.Cin.Pstand) || (Gin = Gin * C.Cin.Pstand)
            return temp * Gin
        elseif ncw == C.k
            base_ring(w) == C.F || throw(ArgumentError("Vector must have the same base ring as the code."))
            return w * C.G
        else
            throw(ArgumentError("Vector has incorrect dimension"))
        end
    else
        len = length(v)
        if len == C.Cout.k
            w = matrix(C.Cout.F, 1, len, v)
            Gout = generatormatrix(C.Cout, true)
            ismissing(C.Cout.Pstand) || (Gout = Gout * C.Cout.Pstand)
            temp = w * Gout
            if C.type == :expanded
                D = _expansiondict(C.Cout.F, C.Cin.F, C.dualbasis)
                temp = _expandmatrix(temp, D, div(degree(C.Cout.F), degree(C.Cin.F)))
                # TODO: this is no longer a vector, only want 1 * temp row instead of full basis?
                temp = temp[1, :]
            end
            # should automatically now be in field of inner code
            Gin = generatormatrix(C.Cin, true)
            ismissing(C.Cin.Pstand) || (Gin = Gin * C.Cin.Pstand)
            return temp * Gin
        elseif len == C.k
            return matrix(C.F, 1, len, v) * C.G
        else
            throw(ArgumentError("Vector has incorrect dimension"))
        end
    end
end

# permute?

#############################
         # Quantum
#############################
