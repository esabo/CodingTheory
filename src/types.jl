# Copyright (c) 2023 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # Classical
#############################

#############################
      # abstract types
#############################

abstract type AbstractCode end
abstract type AbstractNonadditiveCode <: AbstractCode end
abstract type AbstractNonlinearCode <: AbstractCode end
abstract type AbstractAdditiveCode <: AbstractCode end
abstract type AbstractLinearCode <: AbstractAdditiveCode end
abstract type AbstractLDPCCode <: AbstractLinearCode end
abstract type AbstractMatrixProductCode <: AbstractLinearCode end
abstract type AbstractReedMullerCode <: AbstractLinearCode end
abstract type AbstractCyclicCode <: AbstractLinearCode end
abstract type AbstractBCHCode <: AbstractCyclicCode end
abstract type AbstractReedSolomonCode <: AbstractBCHCode end
abstract type AbstractQuasiCyclicCode <: AbstractLinearCode end
abstract type AbstractGeneralizedReedSolomonCode <: AbstractLinearCode end
abstract type AbstractAlgebraicGeometryCode <: AbstractLinearCode end

# const CTMatrixTypes = Union{fpMatrix, fqPolyRepMatrix}
const CTMatrixTypes = MatElem{<:FinFieldElem}

#############################
      # concrete types
#############################

#############################
       # linearcode.jl
#############################

struct WeightEnumerator
    polynomial::Union{fmpz_mpoly, AbstractAlgebra.Generic.MPoly{nf_elem}}
    type::Symbol
end
  
mutable struct LinearCode <: AbstractLinearCode
    F::FqNmodFiniteField # base field
    n::Int # length
    k::Int # dimension
    d::Union{Int, Missing} # minimum distance
    lbound::Int # lower bound on d
    ubound::Int # upper bound on d
    G::fq_nmod_mat
    H::fq_nmod_mat
    Gstand::fq_nmod_mat
    Hstand::fq_nmod_mat
    P::Union{fq_nmod_mat, Missing} # permutation matrix for G -> Gstand
    weightenum::Union{WeightEnumerator, Missing}
end

#############################
         # LDPC.jl
#############################

# TODO: don't like having this here as a subobject - rethink
mutable struct LDPCCode <: AbstractLDPCCode
    C::AbstractLinearCode
    numedges::Int
    vardegs::Vector{Int}
    checkdegs::Vector{Int}
    colbound::Int
    rowbound::Int
    limited::Int
    density::Float64
    isreg::Bool
    tangr::Union{Figure, Missing}
    λ::fmpq_poly
    ρ::fmpq_poly
end

#############################
    # MatrixProductCode.jl
#############################

mutable struct MatrixProductCode <: AbstractMatrixProductCode
    F::FqNmodFiniteField # base field
    n::Int # length
    k::Int # dimension
    d::Union{Int, Missing} # minimum distance
    lbound::Int # lower bound on d
    ubound::Int # upper bound on d
    G::fq_nmod_mat
    H::fq_nmod_mat
    Gstand::fq_nmod_mat
    Hstand::fq_nmod_mat
    P::Union{fq_nmod_mat, Missing} # permutation matrix for G -> Gstand
    weightenum::Union{WeightEnumerator, Missing}
    C::Vector{AbstractLinearCode}
    A::fq_nmod_mat
end

#############################
      # ReedMuller.jl
#############################

mutable struct ReedMullerCode <: AbstractReedMullerCode
    F::Union{FqNmodFiniteField, AbstractAlgebra.GFField{Int64}}
    n::Int # length
    k::Int # dimension
    d::Union{Int, Missing} # minimum distance
    lbound::Int # lower bound on d
    ubound::Int # upper bound on d
    r::Integer # order
    m::Integer # number of variables
    G::Union{gfp_mat, fq_nmod_mat}
    H::Union{gfp_mat, fq_nmod_mat}
    Gstand::Union{gfp_mat, fq_nmod_mat}
    Hstand::Union{gfp_mat, fq_nmod_mat}
    P::Union{fq_nmod_mat, Missing} # permutation matrix for G -> Gstand
    weightenum::Union{WeightEnumerator, Missing}
end

#############################
      # cycliccode.jl
#############################

mutable struct CyclicCode <: AbstractCyclicCode
    F::FqNmodFiniteField # base field
    E::FqNmodFiniteField # splitting field
    R::FqNmodPolyRing # polynomial ring of generator polynomial
    β::fq_nmod # n-th root of primitive element of splitting field
    n::Int # length
    k::Int # dimension
    d::Union{Int, Missing} # minimum distance
    b::Int # offset
    δ::Int # BCH bound
    HT::Int # Hartmann-Tzeng refinement
    lbound::Int # lower bound on d
    ubound::Int # upper bound on d
    qcosets::Vector{Vector{Int}}
    qcosetsreps::Vector{Int}
    defset::Vector{Int}
    g::fq_nmod_poly
    h::fq_nmod_poly
    e::fq_nmod_poly
    G::fq_nmod_mat
    H::fq_nmod_mat
    Gstand::fq_nmod_mat
    Hstand::fq_nmod_mat
    P::Union{fq_nmod_mat, Missing} # permutation matrix for G -> Gstand
    weightenum::Union{WeightEnumerator, Missing}
end
  
mutable struct BCHCode <: AbstractBCHCode
    F::FqNmodFiniteField # base field
    E::FqNmodFiniteField # splitting field
    R::FqNmodPolyRing # polynomial ring of generator polynomial
    β::fq_nmod # n-th root of primitive element of splitting field
    n::Int # length
    k::Int # dimension
    d::Union{Int, Missing} # minimum distance
    b::Int # offset
    δ::Int # BCH bound
    HT::Int # Hartmann-Tzeng refinement
    lbound::Int # lower bound on d
    ubound::Int # upper bound on d
    qcosets::Vector{Vector{Int}}
    qcosetsreps::Vector{Int}
    defset::Vector{Int}
    g::fq_nmod_poly
    h::fq_nmod_poly
    e::fq_nmod_poly
    G::fq_nmod_mat
    H::fq_nmod_mat
    Gstand::fq_nmod_mat
    Hstand::fq_nmod_mat
    P::Union{fq_nmod_mat, Missing} # permutation matrix for G -> Gstand
    weightenum::Union{WeightEnumerator, Missing}
end
  
mutable struct ReedSolomonCode <: AbstractReedSolomonCode
    F::FqNmodFiniteField # base field
    E::FqNmodFiniteField # splitting field
    R::FqNmodPolyRing # polynomial ring of generator polynomial
    β::fq_nmod # n-th root of primitive element of splitting field
    n::Int # length
    k::Int # dimension
    d::Union{Int, Missing} # minimum distance
    b::Int # offset
    δ::Int # BCH bound
    HT::Int # Hartmann-Tzeng refinement
    lbound::Int # lower bound on d
    ubound::Int # upper bound on d
    qcosets::Vector{Vector{Int}}
    qcosetsreps::Vector{Int}
    defset::Vector{Int}
    g::fq_nmod_poly
    h::fq_nmod_poly
    e::fq_nmod_poly
    G::fq_nmod_mat
    H::fq_nmod_mat
    Gstand::fq_nmod_mat
    Hstand::fq_nmod_mat
    P::Union{fq_nmod_mat, Missing} # permutation matrix for G -> Gstand
    weightenum::Union{WeightEnumerator, Missing}
end

#############################
    # quasicycliccode.jl
#############################

mutable struct QuasiCyclicCode <: AbstractQuasiCyclicCode
    F::FqNmodFiniteField # base field
    R::AbstractAlgebra.Generic.ResRing{fq_nmod_poly}
    n::Int # length
    k::Int # dimension
    d::Union{Int, Missing} # minimum distance
    lbound::Int # lower bound on d
    ubound::Int # upper bound on d
    G::Union{fq_nmod_mat, Missing}
    H::Union{fq_nmod_mat, Missing}
    Gstand::Union{fq_nmod_mat, Missing}
    Hstand::Union{fq_nmod_mat, Missing}
    P::Union{fq_nmod_mat, Missing} # permutation matrix for G -> Gstand
    weightenum::Union{WeightEnumerator, Missing}
    l::Int
    m::Int
    A::AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Res{fq_nmod_poly}}
    Atype::Char
    W::Matrix{Int}
    type::Int
end

#############################
 # GeneralizedReedSolomon.jl
#############################

mutable struct GeneralizedReedSolomonCode <: AbstractGeneralizedReedSolomonCode
    F::FqNmodFiniteField # base field
    n::Int # length
    k::Int # dimension
    d::Union{Int, Missing} # minimum distance
    lbound::Int # lower bound on d
    ubound::Int # upper bound on d
    scalars::Vector{fq_nmod}
    dualscalars::Vector{fq_nmod}
    evalpts::Vector{fq_nmod}
    G::fq_nmod_mat
    H::fq_nmod_mat
    Gstand::fq_nmod_mat
    Hstand::fq_nmod_mat
    P::Union{fq_nmod_mat, Missing} # permutation matrix for G -> Gstand
    weightenum::Union{WeightEnumerator, Missing} # TODO: should never be missing? is complete known for MDS?
end

#############################
         # Quantum
#############################

#############################
      # abstract types
#############################

abstract type AbstractSubsystemCode <: AbstractAdditiveCode end
abstract type AbstractSubsystemCodeCSS <: AbstractSubsystemCode end
abstract type AbstractStabilizerCode <: AbstractSubsystemCode end
abstract type AbstractStabilizerCodeCSS <: AbstractStabilizerCode end
abstract type AbstractGraphStateSubsystem <: AbstractSubsystemCode end
abstract type AbstractGraphStateSubsystemCSS <: AbstractSubsystemCodeCSS end
abstract type AbstractGraphStateStabilizer <: AbstractStabilizerCode end
abstract type AbstractGraphStateStabilizerCSS <: AbstractStabilizerCodeCSS end
abstract type AbstractHypergraphProductCode <: AbstractStabilizerCodeCSS end
abstract type AbstractEASubsystemCode <: AbstractSubsystemCode end
abstract type AbstractEASubsystemCodeCSS <: AbstractEASubsystemCode end
abstract type AbstractEAStabilizerCode <: AbstractStabilizerCode end
abstract type AbstractEAStabilizerCodeCSS <: AbstractEAStabilizerCode end

# AbstractQuantumLDPCCode, AbstractQuantumLDPCCSSCode?

#############################
      # concrete types
#############################

#############################
      # subsystemcode.jl
#############################

mutable struct SubsystemCodeCSS <: AbstractSubsystemCodeCSS
    F::FqNmodFiniteField
    n::Int
    k::Union{Int, Rational{BigInt}}
    r::Int
    d::Union{Int, Missing}
    stabs::fq_nmod_mat
    Xstabs::fq_nmod_mat
    Zstabs::fq_nmod_mat
    Xorigcode::Union{LinearCode, Missing}
    ZorigCode::Union{LinearCode, Missing}
    signs::Vector{nmod}
    Xsigns::Vector{nmod}
    Zsigns::Vector{nmod}
    logicals::Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}
    logsmat::fq_nmod_mat
    charvec::Vector{nmod}
    gaugeops::Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}
    gopsmat::fq_nmod_mat
    overcomplete::Bool
end
  
mutable struct SubsystemCode <: AbstractSubsystemCode
    F::FqNmodFiniteField
    n::Int
    k::Union{Int, Rational{BigInt}}
    r::Int
    d::Union{Int, Missing}
    stabs::fq_nmod_mat
    logicals::Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}
    logsmat::fq_nmod_mat
    charvec::Vector{nmod}
    signs::Vector{nmod}
    gaugeops::Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}
    gopsmat::fq_nmod_mat
    overcomplete::Bool
end

#############################
      # stabilizercode.jl
#############################

mutable struct StabilizerCodeCSS <: AbstractStabilizerCodeCSS
    F::FqNmodFiniteField
    n::Int
    k::Union{Int, Rational{BigInt}}
    d::Union{Int, Missing}
    dx::Union{Int, Missing}
    dz::Union{Int, Missing}
    stabs::fq_nmod_mat
    Xstabs::fq_nmod_mat
    Zstabs::fq_nmod_mat
    Xorigcode::Union{LinearCode, Missing}
    ZorigCode::Union{LinearCode, Missing}
    signs::Vector{nmod}
    Xsigns::Vector{nmod}
    Zsigns::Vector{nmod}
    logicals::Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}
    logsmat::fq_nmod_mat
    charvec::Vector{nmod}
    sCWEstabs::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    sCWEdual::Union{WeightEnumerator, Missing} # S^⟂
    sCWElogs::Union{WeightEnumerator, Missing}
    overcomplete::Bool
    pure::Union{Bool, Missing}
end
  
mutable struct StabilizerCode <: AbstractStabilizerCode
    F::FqNmodFiniteField
    n::Int
    k::Union{Int, Rational{BigInt}}
    d::Union{Int, Missing}
    stabs::fq_nmod_mat
    logicals::Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}
    logsmat::fq_nmod_mat
    charvec::Vector{nmod}
    signs::Vector{nmod}
    sCWEstabs::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    sCWEdual::Union{WeightEnumerator, Missing} # S^⟂
    sCWElogs::Union{WeightEnumerator, Missing}
    overcomplete::Bool
    pure::Union{Bool, Missing}
end

#############################
      # graphstate.jl
#############################

mutable struct GraphStateSubsystem <: AbstractGraphStateSubsystem
    F::FqNmodFiniteField
    n::Int
    k::Int
    r::Int
    d::Union{Int, Missing}
    stabs::fq_nmod_mat
    charvec::Vector{nmod}
    signs::Vector{nmod}
    wtenum::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    overcomplete::Bool
    gaugeops::Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}
    gopsmat::fq_nmod_mat
end
  
mutable struct GraphStateSubsystemCSS <: AbstractGraphStateSubsystemCSS
    F::FqNmodFiniteField
    n::Int
    k::Int
    r::Int
    d::Union{Int, Missing}
    dx::Union{Int, Missing}
    dz::Union{Int, Missing}
    stabs::fq_nmod_mat
    Xstabs::fq_nmod_mat
    Zstabs::fq_nmod_mat
    Xorigcode::Union{LinearCode, Missing}
    ZorigCode::Union{LinearCode, Missing}
    signs::Vector{nmod}
    Xsigns::Vector{nmod}
    Zsigns::Vector{nmod}
    charvec::Vector{nmod}
    wtenum::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    overcomplete::Bool
    gaugeops::Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}
    gopsmat::fq_nmod_mat
end
  
mutable struct GraphStateStabilizer <: AbstractGraphStateStabilizer
    F::FqNmodFiniteField
    n::Int
    k::Int
    d::Union{Int, Missing}
    stabs::fq_nmod_mat
    charvec::Vector{nmod}
    signs::Vector{nmod}
    wtenum::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    overcomplete::Bool
end
  
mutable struct GraphStateStabilizerCSS <: AbstractGraphStateStabilizerCSS
    F::FqNmodFiniteField
    n::Int
    k::Int
    d::Union{Int, Missing}
    dx::Union{Int, Missing}
    dz::Union{Int, Missing}
    stabs::fq_nmod_mat
    Xstabs::fq_nmod_mat
    Zstabs::fq_nmod_mat
    Xorigcode::Union{LinearCode, Missing}
    ZorigCode::Union{LinearCode, Missing}
    signs::Vector{nmod}
    Xsigns::Vector{nmod}
    Zsigns::Vector{nmod}
    charvec::Vector{nmod}
    wtenum::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    overcomplete::Bool
end

#############################
   # quantumproductcodes.jl
#############################

# J. Tillich, G. Zémor. "Quantum LDPC codes with positive rate and minimum distance
# proportional to n^(1/2)". (2013) arXiv:0903.0566v2
mutable struct HypergraphProductCode <: AbstractHypergraphProductCode
    F::FqNmodFiniteField
    n::Integer
    k::Union{Integer, Rational{BigInt}}
    d::Union{Integer, Missing}
    dx::Union{Integer, Missing}
    dz::Union{Integer, Missing}
    stabs::fq_nmod_mat
    Xstabs::fq_nmod_mat
    Zstabs::fq_nmod_mat
    C1::Union{LinearCode, Missing}
    C2::Union{LinearCode, Missing}
    signs::Vector{nmod}
    Xsigns::Vector{nmod}
    Zsigns::Vector{nmod}
    dualgens::fq_nmod_mat
    logspace::Union{fq_nmod_mat, Missing}
    logicals::Union{Vector{Tuple{fq_nmod_mat, fq_nmod_mat}}, Missing}
    charvec::Vector{nmod}
    sCWEstabs::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    sCWEdual::Union{WeightEnumerator, Missing} # S^⟂
    overcomplete::Bool
    Lsigns::Union{Vector{nmod}, Missing}
    # TODO: remove dualgens, logspace, missing on logs, Lsigns
end

#############################
         # traits
#############################

const CSSTypes = Union{AbstractSubsystemCodeCSS, AbstractStabilizerCodeCSS, AbstractGraphStateStabilizerCSS, AbstractGraphStateSubsystemCSS}
const GraphStateTypes = Union{AbstractGraphStateSubsystem, AbstractGraphStateSubsystemCSS, AbstractGraphStateStabilizer, AbstractGraphStateStabilizerCSS}

# TODO: base traits on conceptualparent, issubtype
abstract type LogicalTrait end
struct HasLogicals <: LogicalTrait end
struct HasNoLogicals <: LogicalTrait end
LogicalTrait(::Type{T}) where {T <: AbstractSubsystemCode} = HasLogicals()
LogicalTrait(::Type{T}) where {T <: GraphStateTypes} = HasNoLogicals()

abstract type GaugeTrait end
struct HasGauges <: GaugeTrait end
struct HasNoGauges <: GaugeTrait end
GaugeTrait(::Type{T}) where {T <: AbstractSubsystemCode} = HasGauges()
GaugeTrait(::Type{T}) where {T <: AbstractStabilizerCode} = HasNoGauges()

abstract type CSSTrait end
struct IsCSS <: CSSTrait end
struct IsNotCSS <: CSSTrait end
CSSTrait(::Type{T}) where {T <: AbstractSubsystemCode} = IsNotCSS()
CSSTrait(::Type{T}) where {T <: CSSTypes} = IsCSS()
