# Copyright (c) 2025 Eric Sabo, David Marquis
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

# TODO do we even want to accept a Gram matrix here?
"""
    Lattice(B::MatrixElem{<:RingElem}; type::Symbol = :basis)

Return a lattice using the columns of `B` as a basis.
"""
function Lattice(B::MatrixElem{<:RingElem}; type::Symbol = :basis)
    if type == :basis
        rank(B) == ncols(B) || throw(ArgumentError("Basis must have full column rank."))
        return Lattice(B)
    elseif type == :Gram
        # which decomp to use here?

        return Lattice(B)
    else
        throw(ArgumentError("Parameter `type` must be either `:basis` or `:Gram`."))
    end
end

"""
    Lattice(V::Vector{MatrixElem{<:RingElem}})

Return a lattice using the elements of the vectors as a basis.
"""
function Lattice(V::Vector{MatrixElem{<:RingElem}})
    # check equal sizes and base rings
    # stack rows into columns of matrix

    rank(B) == ncols(B) || throw(ArgumentError("Basis must have full column rank."))
    return Lattice(B)
end

"""


"""
function Lattice(G::MatrixElem{<:RingElem})


    return Lattice(B)
end

"""
    construction_A(C::AbstractLinearCode)

Return the lattice given by applying Construction A to the code `C`.
"""
function construction_A(C::AbstractLinearCode)


end

"""
    construction_B(C::AbstractLinearCode)

Return the lattice given by applying Construction B to the code `C`.
"""
function construction_B(C::AbstractLinearCode)


end

"""
    construction_C(C::AbstractLinearCode)

Return the lattice given by applying Construction C to the code `C`.
"""
function construction_C(C::AbstractLinearCode)


end

#############################
      # getter functions
#############################

"""
    basis(L::Lattice)

Return the basis of the lattice.
"""
basis(L::Lattice) = L.B

"""
    rank(L::Lattice)

Return the rank of the lattice.
"""
rank(L::Lattice) = size(L.B, 1)

"""
    dimension(L::Lattice)
    dim(L::Lattice)

Return the dimension of the lattice.
"""
dimension(L::Lattice) = size(L.B, 2)
dim(L::Lattice) = dimension(L)

"""
    is_full_rank(L::Lattice)

Return `true` if the lattice is full rank; otherwise return `false`.
"""
is_full_rank(L::Lattice) = dimension(L) == rank(L)

"""
    discriminant(L::Lattice)
    disc(L::Lattice)

Return the discriminant of the lattice.
"""
discriminant(L::Lattice) = det(Gram_matrix(L))
disc(L::Lattice) = discriminant(L)

"""
    Gram_matrix(L::Lattice)

Return the Gram matrix of the basis of the lattice.
"""
Gram_matrix(L::Lattice) = gram(L.B)

"""
    volume(L::Lattice)
    vol(L::Lattice)
    determinate(L::Lattice)
    det(L::Lattice)

Return the volumne of the lattice.
"""
volume(L::Lattice) = sqrt(disc(L))
vol(L::Lattice) = volume(L)
determinate(L::Lattice) = volume(L)
det(L::Lattice) = volume(L)

"""
    check_matrix(L::Lattice)

Return the check matrix of the lattice.
"""
check_matrix(L::Lattice) = inverse(L.B)

"""
    shortest_vector(L::Lattice)

"""
shortest_vector(L::Lattice)

"""
    successive_minima(L::Lattice)

"""
successive_minima(L::Lattice)

"""
    dual_basis(L::Lattice)

Return a dual basis for the lattice.
"""
function dual_basis(L::Lattice)
    Ginv = inverse(Gram_matrix(L))
    new_basis = Ginv * B
    # b^*_i = \sum_j (G_inv)_ij b_j
    return new_basis
end

"""
    dual(L::Lattice)

Return the dual lattice.
"""
dual(L::Lattice) = Lattice(dual_basis(L))

"""
    is_unimodular(L::Lattice)

Return `true` if the lattie is unimodular; otherwise return `false`.
"""
function is_unimodular(L::Lattice)
    # Gram matrix has integer entries and disc(L) == 1
end

#############################
     # general functions
#############################

"""


"""
function in(v::MatrixElem{<:RingElem}, L::Lattice)
    # check sizes and rings
    # setup Bx = v
    # true if solution exists and is in ZZ
end

function is_sublattice(L1::Lattice, L2::Lattice)

end

function is_equal(L1::Lattice, L2::Lattice)
    # L1 ⊆ L2 && L2 ⊆ L1

end

"""


"""
function coding_gain(L::Lattice)
    
end

"""


"""
function covering_radius(L::Lattice)

end

# SVP, CVP
