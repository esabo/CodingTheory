# Copyright (c) 2023 Eric Sabo, Benjamin Ide
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

# easy to change later if necessary:
const CTFieldTypes = FinField
const CTFieldElem = FinFieldElem
const CTMatrixTypes = MatElem{<:CTFieldElem}
const CTPolyRing = PolyRing{<:CTFieldElem}
const CTPolyRingElem = PolyRingElem{<:CTFieldElem}

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
    F::CTFieldTypes # base field
    n::Int # length
    k::Int # dimension
    d::Union{Int, Missing} # minimum distance
    lbound::Int # lower bound on d
    ubound::Int # upper bound on d
    G::CTMatrixTypes
    H::CTMatrixTypes
    Gstand::CTMatrixTypes
    Hstand::CTMatrixTypes
    Pstand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> Gstand
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
    F::CTFieldTypes # base field
    n::Int # length
    k::Int # dimension
    d::Union{Int, Missing} # minimum distance
    lbound::Int # lower bound on d
    ubound::Int # upper bound on d
    G::CTMatrixTypes
    H::CTMatrixTypes
    Gstand::CTMatrixTypes
    Hstand::CTMatrixTypes
    Pstand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> Gstand
    weightenum::Union{WeightEnumerator, Missing}
    C::Vector{AbstractLinearCode}
    A::fq_nmod_mat
end

#############################
      # ReedMuller.jl
#############################

mutable struct ReedMullerCode <: AbstractReedMullerCode
    F::CTFieldTypes
    n::Int # length
    k::Int # dimension
    d::Union{Int, Missing} # minimum distance
    lbound::Int # lower bound on d
    ubound::Int # upper bound on d
    r::Integer # order
    m::Integer # number of variables
    G::CTMatrixTypes
    H::CTMatrixTypes
    Gstand::CTMatrixTypes
    Hstand::CTMatrixTypes
    Pstand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> Gstand
    weightenum::Union{WeightEnumerator, Missing}
end

#############################
      # cycliccode.jl
#############################

mutable struct CyclicCode <: AbstractCyclicCode
    F::CTFieldTypes # base field
    E::CTFieldTypes # splitting field
    R::CTPolyRing # polynomial ring of generator polynomial
    β::CTFieldElem # n-th root of primitive element of splitting field
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
    g::CTPolyRingElem
    h::CTPolyRingElem
    e::CTPolyRingElem
    G::CTMatrixTypes
    H::CTMatrixTypes
    Gstand::CTMatrixTypes
    Hstand::CTMatrixTypes
    Pstand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> Gstand
    weightenum::Union{WeightEnumerator, Missing}
end
  
mutable struct BCHCode <: AbstractBCHCode
    F::CTFieldTypes # base field
    E::CTFieldTypes # splitting field
    R::CTPolyRing # polynomial ring of generator polynomial
    β::CTFieldElem # n-th root of primitive element of splitting field
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
    g::CTPolyRingElem
    h::CTPolyRingElem
    e::CTPolyRingElem
    G::CTMatrixTypes
    H::CTMatrixTypes
    Gstand::CTMatrixTypes
    Hstand::CTMatrixTypes
    Pstand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> Gstand
    weightenum::Union{WeightEnumerator, Missing}
end
  
mutable struct ReedSolomonCode <: AbstractReedSolomonCode
    F::CTFieldTypes # base field
    E::CTFieldTypes # splitting field
    R::CTPolyRing # polynomial ring of generator polynomial
    β::CTFieldElem # n-th root of primitive element of splitting field
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
    g::CTPolyRingElem
    h::CTPolyRingElem
    e::CTPolyRingElem
    G::CTMatrixTypes
    H::CTMatrixTypes
    Gstand::CTMatrixTypes
    Hstand::CTMatrixTypes
    Pstand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> Gstand
    weightenum::Union{WeightEnumerator, Missing}
end

#############################
    # quasicycliccode.jl
#############################

mutable struct QuasiCyclicCode <: AbstractQuasiCyclicCode
    F::CTFieldTypes # base field
    R::AbstractAlgebra.Generic.ResRing
    n::Int # length
    k::Int # dimension
    d::Union{Int, Missing} # minimum distance
    lbound::Int # lower bound on d
    ubound::Int # upper bound on d
    G::Union{CTMatrixTypes, Missing}
    H::Union{CTMatrixTypes, Missing}
    Gstand::Union{CTMatrixTypes, Missing}
    Hstand::Union{CTMatrixTypes, Missing}
    Pstand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> Gstand
    weightenum::Union{WeightEnumerator, Missing}
    l::Int
    m::Int
    A::MatElem{<:ResElem}
    Atype::Symbol
    W::Matrix{Int}
    type::Int
end

#############################
 # GeneralizedReedSolomon.jl
#############################

mutable struct GeneralizedReedSolomonCode <: AbstractGeneralizedReedSolomonCode
    F::CTFieldTypes # base field
    n::Int # length
    k::Int # dimension
    d::Union{Int, Missing} # minimum distance
    lbound::Int # lower bound on d
    ubound::Int # upper bound on d
    scalars::Vector{<:CTFieldElem}
    dualscalars::Vector{<:CTFieldElem}
    evalpts::Vector{<:CTFieldElem}
    G::CTMatrixTypes
    H::CTMatrixTypes
    Gstand::CTMatrixTypes
    Hstand::CTMatrixTypes
    Pstand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> Gstand
    weightenum::Union{WeightEnumerator, Missing} # TODO: should never be missing? is completely known for MDS?
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
    F::CTFieldTypes
    n::Int
    k::Union{Int, Rational{BigInt}}
    r::Int
    d::Union{Int, Missing}
    stabs::CTMatrixTypes
    Xstabs::CTMatrixTypes
    Zstabs::CTMatrixTypes
    Xorigcode::Union{LinearCode, Missing}
    ZorigCode::Union{LinearCode, Missing}
    signs::Vector{nmod}
    Xsigns::Vector{nmod}
    Zsigns::Vector{nmod}
    logicals::Vector{Tuple{T, T}} where T <: CTMatrixTypes
    logsmat::CTMatrixTypes
    charvec::Vector{nmod}
    gaugeops::Vector{Tuple{T, T}} where T <: CTMatrixTypes
    gopsmat::CTMatrixTypes
    overcomplete::Bool
    stabsstand::CTMatrixTypes
    standr::Int
    standk::Int
    Pstand::Union{CTMatrixTypes, Missing}
end
  
mutable struct SubsystemCode <: AbstractSubsystemCode
    F::CTFieldTypes
    n::Int
    k::Union{Int, Rational{BigInt}}
    r::Int
    d::Union{Int, Missing}
    stabs::CTMatrixTypes
    logicals::Vector{Tuple{T, T}} where T <: CTMatrixTypes
    logsmat::CTMatrixTypes
    charvec::Vector{nmod}
    signs::Vector{nmod}
    gaugeops::Vector{Tuple{T, T}} where T <: CTMatrixTypes
    gopsmat::CTMatrixTypes
    overcomplete::Bool
    stabsstand::CTMatrixTypes
    standr::Int
    standk::Int
    Pstand::Union{CTMatrixTypes, Missing}
end

#############################
      # stabilizercode.jl
#############################

mutable struct StabilizerCodeCSS <: AbstractStabilizerCodeCSS
    F::CTFieldTypes
    n::Int
    k::Union{Int, Rational{BigInt}}
    d::Union{Int, Missing}
    dx::Union{Int, Missing}
    dz::Union{Int, Missing}
    stabs::CTMatrixTypes
    Xstabs::CTMatrixTypes
    Zstabs::CTMatrixTypes
    Xorigcode::Union{LinearCode, Missing}
    ZorigCode::Union{LinearCode, Missing}
    signs::Vector{nmod}
    Xsigns::Vector{nmod}
    Zsigns::Vector{nmod}
    logicals::Vector{Tuple{T, T}} where T <: CTMatrixTypes
    logsmat::CTMatrixTypes
    charvec::Vector{nmod}
    sCWEstabs::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    sCWEdual::Union{WeightEnumerator, Missing} # S^⟂
    sCWElogs::Union{WeightEnumerator, Missing}
    overcomplete::Bool
    pure::Union{Bool, Missing}
    stabsstand::CTMatrixTypes
    standr::Int
    standk::Int
    Pstand::Union{CTMatrixTypes, Missing}
end
  
mutable struct StabilizerCode <: AbstractStabilizerCode
    F::CTFieldTypes
    n::Int
    k::Union{Int, Rational{BigInt}}
    d::Union{Int, Missing}
    stabs::CTMatrixTypes
    logicals::Vector{Tuple{T, T}} where T <: CTMatrixTypes
    logsmat::CTMatrixTypes
    charvec::Vector{nmod}
    signs::Vector{nmod}
    sCWEstabs::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    sCWEdual::Union{WeightEnumerator, Missing} # S^⟂
    sCWElogs::Union{WeightEnumerator, Missing}
    overcomplete::Bool
    pure::Union{Bool, Missing}
    stabsstand::CTMatrixTypes
    standr::Int
    standk::Int
    Pstand::Union{CTMatrixTypes, Missing}
end

#############################
      # graphstate.jl
#############################

mutable struct GraphStateSubsystem <: AbstractGraphStateSubsystem
    F::CTFieldTypes
    n::Int
    k::Int
    r::Int
    d::Union{Int, Missing}
    stabs::CTMatrixTypes
    charvec::Vector{nmod}
    signs::Vector{nmod}
    wtenum::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    overcomplete::Bool
    gaugeops::Vector{Tuple{T, T}} where T <: CTMatrixTypes
    gopsmat::CTMatrixTypes
    stabsstand::CTMatrixTypes
    standr::Int
    standk::Int
    Pstand::Union{CTMatrixTypes, Missing}
end
  
mutable struct GraphStateSubsystemCSS <: AbstractGraphStateSubsystemCSS
    F::CTFieldTypes
    n::Int
    k::Int
    r::Int
    d::Union{Int, Missing}
    dx::Union{Int, Missing}
    dz::Union{Int, Missing}
    stabs::CTMatrixTypes
    Xstabs::CTMatrixTypes
    Zstabs::CTMatrixTypes
    Xorigcode::Union{LinearCode, Missing}
    ZorigCode::Union{LinearCode, Missing}
    signs::Vector{nmod}
    Xsigns::Vector{nmod}
    Zsigns::Vector{nmod}
    charvec::Vector{nmod}
    wtenum::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    overcomplete::Bool
    gaugeops::Vector{Tuple{T, T}} where T <: CTMatrixTypes
    gopsmat::CTMatrixTypes
    stabsstand::CTMatrixTypes
    standr::Int
    standk::Int
    Pstand::Union{CTMatrixTypes, Missing}
end
  
mutable struct GraphStateStabilizer <: AbstractGraphStateStabilizer
    F::CTFieldTypes
    n::Int
    k::Int
    d::Union{Int, Missing}
    stabs::CTMatrixTypes
    charvec::Vector{nmod}
    signs::Vector{nmod}
    wtenum::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    overcomplete::Bool
    stabsstand::CTMatrixTypes
    standr::Int
    standk::Int
    Pstand::Union{CTMatrixTypes, Missing}
end
  
mutable struct GraphStateStabilizerCSS <: AbstractGraphStateStabilizerCSS
    F::CTFieldTypes
    n::Int
    k::Int
    d::Union{Int, Missing}
    dx::Union{Int, Missing}
    dz::Union{Int, Missing}
    stabs::CTMatrixTypes
    Xstabs::CTMatrixTypes
    Zstabs::CTMatrixTypes
    Xorigcode::Union{LinearCode, Missing}
    ZorigCode::Union{LinearCode, Missing}
    signs::Vector{nmod}
    Xsigns::Vector{nmod}
    Zsigns::Vector{nmod}
    charvec::Vector{nmod}
    wtenum::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    overcomplete::Bool
    stabsstand::CTMatrixTypes
    standr::Int
    standk::Int
    Pstand::Union{CTMatrixTypes, Missing}
end

#############################
   # quantumproductcodes.jl
#############################

# J. Tillich, G. Zémor. "Quantum LDPC codes with positive rate and minimum distance
# proportional to n^(1/2)". (2013) arXiv:0903.0566v2
mutable struct HypergraphProductCode <: AbstractHypergraphProductCode
    F::CTFieldTypes
    n::Integer
    k::Union{Integer, Rational{BigInt}}
    d::Union{Integer, Missing}
    dx::Union{Integer, Missing}
    dz::Union{Integer, Missing}
    stabs::CTMatrixTypes
    Xstabs::CTMatrixTypes
    Zstabs::CTMatrixTypes
    C1::Union{LinearCode, Missing}
    C2::Union{LinearCode, Missing}
    signs::Vector{nmod}
    Xsigns::Vector{nmod}
    Zsigns::Vector{nmod}
    logicals::Vector{Tuple{T, T}} where T <: CTMatrixTypes
    logsmat::CTMatrixTypes
    charvec::Vector{nmod}
    overcomplete::Bool
    stabsstand::CTMatrixTypes
    standr::Int
    standk::Int
    Pstand::Union{CTMatrixTypes, Missing}
    sCWEstabs::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    sCWEdual::Union{WeightEnumerator, Missing} # S^⟂
end

#############################
         # traits
#############################

const CSSTypes = Union{AbstractSubsystemCodeCSS, AbstractStabilizerCodeCSS, AbstractGraphStateStabilizerCSS, AbstractGraphStateSubsystemCSS, AbstractHypergraphProductCode}
const GraphStateTypes = Union{AbstractGraphStateSubsystem, AbstractGraphStateSubsystemCSS, AbstractGraphStateStabilizer, AbstractGraphStateStabilizerCSS}

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

#############################
       # copy function
#############################

"""
    copy(C::T) where T <: AbstractCode

Returns a copy of the code `C`.
"""
function copy(C::T) where T <: AbstractCode
    C2 = deepcopy(C)
    hasfield(T, :F) && (C2.F = C.F;)
    hasfield(T, :E) && (C2.E = C.E;)
    hasfield(T, :R) && (C2.R = C.R;)
    if hasfield(T, :C)
        if isa(C, AbstractLDPCCode)
            S = typeof(C.C)
            hasfield(S, :F) && (C2.C.F = C.C.F;)
            hasfield(S, :E) && (C2.C.E = C.C.E;)
            hasfield(S, :R) && (C2.C.R = C.C.R;)
        elseif isa(C, AbstractMatrixProductCode)
            for i in eachindex(C.C)
                S = typeof(C.C[i])
                hasfield(S, :F) && (C2.C[i].F = C.C[i].F;)
                hasfield(S, :E) && (C2.C[i].E = C.C[i].E;)
                hasfield(S, :R) && (C2.C[i].R = C.C[i].R;)
            end
        end
        hasfield(S, :C) && @warn "Some sub-sub-codes in the copied struct will have deepcopied Galois fields"
    end
    return C2
end
