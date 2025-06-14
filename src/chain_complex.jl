# Copyright (c) 2023 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
# constructors
#############################

# Base.zero(phi::AbstractAlgebra.Generic.ModuleHomomorphism) = hom(domain(phi), codomain(phi), [zero(codomain(phi)) for i in 1:ngens(domain(phi))])

"""
    chain_complex(∂s::Vector{T}; lowest_degree::Int=0) where T <: CTMatrixTypes

Return the chain complex whose boundary maps are given by `∂s` and lowest degree by
`lowest_degree`.
"""
function chain_complex(∂s::Vector{T}; lowest_degree::Int = 0) where {T<:CTMatrixTypes}
    F = base_ring(∂s[1])
    all(x -> base_ring(x) == F, ∂s) ||
        throw(ArgumentError("All inputs must be over the same base ring"))
    spaces = [[vector_space(F, ncols(∂i)) for ∂i in ∂s]; vector_space(F, nrows(∂s[end]))]
    # need to transpose because they are forced into right multiplication for noncommuntative rings
    morphisms = [hom(spaces[i], spaces[i+1], transpose(∂i)) for (i, ∂i) in enumerate(∂s)]
    # ComplexOfMorphisms(typeof(spaces[1]), morphisms, typ = :chain, seed = lowest_degree)
    return chain_complex(morphisms; seed = lowest_degree)
end

# these are not defined for subsystem CSS codes or graph states?
"""
    chain_complex(S::AbstractStabilizerCode)

Return the 3-term chain complex associated with the CSS code `S`.
"""
chain_complex(S::AbstractStabilizerCode) = chain_complex(CSSTrait(typeof(S)), S)
chain_complex(::IsCSS, S::AbstractStabilizerCode) =
    chain_complex([transpose(S.Z_stabs), S.X_stabs])
chain_complex(::IsNotCSS, S::AbstractStabilizerCode) =
    throw(ArgumentError("This is only defined for CSS codes"))

"""
    chain_complex(C::AbstractLinearCode)

Return the 2-term chain complex associated with the parity-check matrix of `C`.
"""
chain_complex(C::AbstractLinearCode) = chain_complex([parity_check_matrix(C)])

#############################
# getter functions
#############################

"""
    field(chain::CTChainComplex)

Return the base ring of the boundary maps of the chain complex.
"""
field(chain::CTChainComplex) = base_ring(chain[1])

"""
    length(chain::CTChainComplex)

Return the number of spaces in the chain complex.

* Note
- This is the number of boundary maps plus one.
"""
length(chain::CTChainComplex) = length(chain.maps) + 1

"""
    boundaries(chain::CTChainComplex)

Return the boundary maps of the chain complex.

* Note
- boundaries(chain)[i] is the map \$\\partial_i: C_i \\to C_{i - 1}\$.
- The underlying chain complexes are built around right multiplication; this casts
  the results back to the cmmonly assumed left multiplication maps.
"""
boundaries(chain::CTChainComplex) = map(transpose, chain.maps)

"""
    type(chain::CTChainComplex)

Return `:chain` if `chain` is a chain complex or `:cochain` if it is a cochain complex.
"""
type(chain::CTChainComplex) = chain.typ

#############################
# setter functions
#############################

#############################
# general functions
#############################

# TODO: redo
"""
    cochain(chain::ChainComplex)

Return the dual of the chain complex.
"""
cochain(chain::ChainComplex) =
    ChainComplex([transpose(∂) for ∂ in reverse(chain.boundaries)])

"""
    ⊗(chain_A::ChainComplex{T}, chain_B::ChainComplex{T}) where T <: CTMatrixTypes
    tensor_product(chain_A::ChainComplex{T}, chain_B::ChainComplex{T}) where T <: CTMatrixTypes

Return the total complex of the tensor product of the two chains.
"""
⊗(chain_A::CTChainComplex, chain_B::CTChainComplex) = tensor_product(chain_A, chain_B)
double_complex(chain_A::CTChainComplex, chain_B::CTChainComplex) = chain_A ⊗ chain_B
bicomplex(chain_A::CTChainComplex, chain_B::CTChainComplex) = chain_A ⊗ chain_B


double_complex(C1::AbstractLinearCode, C2::AbstractLinearCode) =
    chain_complex(C1) ⊗ chain_complex(C2)
double_complex(S::AbstractStabilizerCodeCSS, C::AbstractLinearCode) =
    chain_complex(S) ⊗ chain_complex(C)
double_complex(S1::AbstractStabilizerCodeCSS, S2::AbstractStabilizerCodeCSS) =
    chain_complex(S1) ⊗ chain_complex(S2)

total_complex(C1::AbstractLinearCode, C2::AbstractLinearCode) =
    total_complex(double_complex(C1, C2))
total_complex(S::AbstractStabilizerCodeCSS, C::AbstractLinearCode) =
    total_complex(double_complex(S, C))
total_complex(S1::AbstractStabilizerCodeCSS, S2::AbstractStabilizerCodeCSS) =
    total_complex(double_complex(S1, S2))

"""
    distance_balancing(S::StabilizerCodeCSS, C::AbstractLinearCode)

Return the distance-balanced code of `S` and `C`.
"""
function distance_balancing(S::StabilizerCodeCSS, C::AbstractLinearCode)
    is_overcomplete(C, :H) && throw(ArgumentError("Parity-check matrix is overcomplete"))

    ∂ = boundaries(tensor_product(S, C))
    return CSSCode(transpose(∂[1]), ∂[2])
end
