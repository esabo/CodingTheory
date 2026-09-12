# Copyright (c) 2023 - 2024 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
      # abstract types
#############################

"""
$(TYPEDEF)

Supertype for additive quantum subsystem codes, including stabilizer codes as the gauge-free specialization.
"""
abstract type AbstractSubsystemCode <: AbstractAdditiveCode end

"""
$(TYPEDEF)

Supertype for CSS subsystem codes with separate ``X``- and ``Z``-type stabilizer data.
"""
abstract type AbstractSubsystemCodeCSS <: AbstractSubsystemCode end

"""
$(TYPEDEF)

Supertype for stabilizer codes, represented in this hierarchy as subsystem codes without gauge qubits.
"""
abstract type AbstractStabilizerCode <: AbstractSubsystemCode end

"""
$(TYPEDEF)

Supertype for CSS stabilizer codes with separate ``X``- and ``Z``-type stabilizer data.
"""
abstract type AbstractStabilizerCodeCSS <: AbstractStabilizerCode end

"""
$(TYPEDEF)

Supertype for graph-state subsystem-code representations.
"""
abstract type AbstractGraphStateSubsystem <: AbstractSubsystemCode end

"""
$(TYPEDEF)

Supertype for CSS graph-state subsystem-code representations.
"""
abstract type AbstractGraphStateSubsystemCSS <: AbstractSubsystemCodeCSS end

"""
$(TYPEDEF)

Supertype for graph-state stabilizer-code representations.
"""
abstract type AbstractGraphStateStabilizer <: AbstractStabilizerCode end

"""
$(TYPEDEF)

Supertype for CSS graph-state stabilizer-code representations.
"""
abstract type AbstractGraphStateStabilizerCSS <: AbstractStabilizerCodeCSS end

"""
$(TYPEDEF)

Supertype for CSS stabilizer codes obtained from the hypergraph-product construction.
"""
abstract type AbstractHypergraphProductCode <: AbstractStabilizerCodeCSS end

"""
$(TYPEDEF)

Supertype for entanglement-assisted subsystem codes.
"""
abstract type AbstractEASubsystemCode <: AbstractSubsystemCode end

"""
$(TYPEDEF)

Supertype for CSS entanglement-assisted subsystem codes.
"""
abstract type AbstractEASubsystemCodeCSS <: AbstractEASubsystemCode end

"""
$(TYPEDEF)

Supertype for entanglement-assisted stabilizer codes.
"""
abstract type AbstractEAStabilizerCode <: AbstractStabilizerCode end

"""
$(TYPEDEF)

Supertype for CSS entanglement-assisted stabilizer codes.
"""
abstract type AbstractEAStabilizerCodeCSS <: AbstractEAStabilizerCode end
abstract type AbstractMonomialCode <: AbstractStabilizerCodeCSS end
abstract type AbstractBivariateBicycleCode <: AbstractMonomialCode end
abstract type AbstractGeneralizedToricCode <: AbstractMonomialCode end
abstract type AbstractGeneralized3DToricCode <: AbstractMonomialCode end

# AbstractQuantumLDPCCode, AbstractQuantumLDPCCSSCode?

abstract type AbstractQuantumNoiseChannel <: AbstractChannel end

function Base.getproperty(S::AbstractSubsystemCode, property::Symbol)
    property in fieldnames(typeof(S)) &&
        return getfield(S, property)
    cache = getfield(S, :cache)
    haskey(cache, property) && return cache[property]
    property == :weight_enum && return weight_enumerator(S)
    property == :weight_dist && return weight_distribution(S)
    throw(ErrorException(
        "type $(typeof(S)) has no field or cached property `$property`."))
end

function Base.setproperty!(
    S::AbstractSubsystemCode, property::Symbol, value
)
    if property in fieldnames(typeof(S))
        return setfield!(S, property, value)
    end
    getfield(S, :cache)[property] = value
    return value
end

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
    X_stabs::CTMatrixTypes
    Z_stabs::CTMatrixTypes
    gauge_ops::Vector{Tuple{CTMatrixTypes, CTMatrixTypes}}
    char_vec::Vector{zzModRingElem}
    cache::Dict{Symbol, Any}
end

mutable struct SubsystemCode <: AbstractSubsystemCode
    F::CTFieldTypes
    n::Int
    k::Union{Int, Rational{BigInt}}
    r::Int
    stabs::CTMatrixTypes
    gauge_ops::Vector{Tuple{CTMatrixTypes, CTMatrixTypes}}
    char_vec::Vector{zzModRingElem}
    cache::Dict{Symbol, Any}
end

#############################
      # stabilizercode.jl
#############################

mutable struct StabilizerCodeCSS <: AbstractStabilizerCodeCSS
    F::CTFieldTypes
    n::Int
    k::Union{Int, Rational{BigInt}}
    X_stabs::CTMatrixTypes
    Z_stabs::CTMatrixTypes
    char_vec::Vector{zzModRingElem}
    cache::Dict{Symbol, Any}
end

mutable struct StabilizerCode <: AbstractStabilizerCode
    F::CTFieldTypes
    n::Int
    k::Union{Int, Rational{BigInt}}
    stabs::CTMatrixTypes
    char_vec::Vector{zzModRingElem}
    cache::Dict{Symbol, Any}
end

#############################
      # graphstate.jl
#############################

mutable struct GraphStateSubsystem <: AbstractGraphStateSubsystem
    F::CTFieldTypes
    n::Int
    k::Int
    r::Int
    stabs::CTMatrixTypes
    gauge_ops::Vector{Tuple{CTMatrixTypes, CTMatrixTypes}}
    char_vec::Vector{zzModRingElem}
    cache::Dict{Symbol, Any}
end

mutable struct GraphStateSubsystemCSS <: AbstractGraphStateSubsystemCSS
    F::CTFieldTypes
    n::Int
    k::Int
    r::Int
    X_stabs::CTMatrixTypes
    Z_stabs::CTMatrixTypes
    gauge_ops::Vector{Tuple{CTMatrixTypes, CTMatrixTypes}}
    char_vec::Vector{zzModRingElem}
    cache::Dict{Symbol, Any}
end

mutable struct GraphStateStabilizer <: AbstractGraphStateStabilizer
    F::CTFieldTypes
    n::Int
    k::Int
    stabs::CTMatrixTypes
    char_vec::Vector{zzModRingElem}
    cache::Dict{Symbol, Any}
end

mutable struct GraphStateStabilizerCSS <: AbstractGraphStateStabilizerCSS
    F::CTFieldTypes
    n::Int
    k::Int
    X_stabs::CTMatrixTypes
    Z_stabs::CTMatrixTypes
    char_vec::Vector{zzModRingElem}
    cache::Dict{Symbol, Any}
end

#############################
# Quantum code-family implementations
#############################


#############################
# Quantum/BB_codes.jl
#############################

"""
$(TYPEDEF)

An algebraic bivariate-bicycle datum before a finite lattice is chosen.

This is deliberately not an `AbstractSubsystemCode`: it has no finite block
length or stabilizer matrix.
"""
struct InfiniteBBCode{T, U, V}
    LR::T
    F::CTFieldTypes
    a::U
    b::V
end

"""
$(TYPEDEF)

A finite bivariate-bicycle CSS stabilizer code.

Representation-specific information (a standard quotient, a twisted Laurent
lattice, or a coprime univariate quotient) is retained in `R`, `a1`, `a2`, and
`N`. Derived matrices and code metadata live in `cache`.
"""
mutable struct BBCode{T, U, V} <: AbstractStabilizerCodeCSS
    R::T
    F::CTFieldTypes
    a::U
    b::V
    a1::Union{Nothing, Tuple{Int, Int}}
    a2::Union{Nothing, Tuple{Int, Int}}
    N::Union{Nothing, Int}
    n::Int
    k::Int
    twisted::Bool
    representation::Symbol
    cache::Dict{Symbol, Any}
end

#############################
# Quantum/generalized_3d_toric_codes.jl
#############################

"""
$(TYPEDEF)

An algebraic three-dimensional generalized toric-code datum.
"""
struct Generalized3DToricCode{T, U, V}
    LR::T
    F::CTFieldTypes
    a::U
    b::V
end

"""
$(TYPEDEF)

A finite three-dimensional generalized toric CSS stabilizer code.
"""
mutable struct FiniteGeneralized3DToricCode{T, U, V} <: AbstractStabilizerCodeCSS
    LR::T
    F::CTFieldTypes
    a::U
    b::V
    a1::Tuple{Int, Int}
    a2::Tuple{Int, Int}
    l_z::Int
    n::Int
    k::Int
    twisted::Bool
    cache::Dict{Symbol, Any}
end

"""
$(TYPEDEF)

A stabilizer code formed by concatenating an outer stabilizer code with an inner stabilizer code.
"""
mutable struct QuantumConcatenatedCode <: AbstractStabilizerCode
    F::CTFieldTypes
    outer_code::AbstractStabilizerCode
    inner_code::AbstractStabilizerCode
    n::Int
    k::Int
    char_vec::Vector{zzModRingElem}
    cache::Dict{Symbol, Any}
end

mutable struct HypergraphProductCode <: AbstractHypergraphProductCode
    F::CTFieldTypes
    C1::AbstractLinearCode
    C2::AbstractLinearCode
    C1T::AbstractLinearCode
    C2T::AbstractLinearCode
    n::Int
    k::Int
    char_vec::Vector{zzModRingElem}
    cache::Dict{Symbol, Any}
end

"""
$(TYPEDEF)

A generalized Shor subsystem code constructed from two classical linear codes.
"""
mutable struct GeneralizedShorCode <: AbstractSubsystemCode
    F::CTFieldTypes
    C1::AbstractLinearCode
    C2::AbstractLinearCode
    n::Int
    k::Int
    r::Int
    char_vec::Vector{zzModRingElem}
    cache::Dict{Symbol, Any}
end

mutable struct HyperBicycleCodeCSS{T <: CTMatrixTypes} <: AbstractStabilizerCodeCSS
    F::CTFieldTypes
    a::Vector{T}
    b::Vector{T}
    χ::Int
    n::Int
    k::Int
    char_vec::Vector{zzModRingElem}
    cache::Dict{Symbol, Any}
end

mutable struct HyperBicycleCode{T <: CTMatrixTypes} <: AbstractStabilizerCode
    F::CTFieldTypes
    a::Vector{T}
    b::Vector{T}
    χ::Int
    n::Int
    k::Int
    char_vec::Vector{zzModRingElem}
    cache::Dict{Symbol, Any}
end

mutable struct GeneralizedBicycleCode{T <: CTMatrixTypes} <: AbstractStabilizerCodeCSS
    F::CTFieldTypes
    A::T
    B::T
    n::Int
    k::Int
    char_vec::Vector{zzModRingElem}
    cache::Dict{Symbol, Any}
end

"""
$(TYPEDEF)

A CSS stabilizer code defined by a lifted-product construction from matrices `A` and `B`.
"""
mutable struct LiftedProductCode{T} <: AbstractStabilizerCodeCSS
    F::CTFieldTypes
    A::T
    B::T
    n::Int
    k::Int
    char_vec::Vector{zzModRingElem}
    cache::Dict{Symbol, Any}
end

"""
$(TYPEDEF)

A bias-tailored lifted-product stabilizer code defined by matrices `A` and `B`.
"""
mutable struct BiasTailoredLiftedProductCode{T} <: AbstractStabilizerCode
    F::CTFieldTypes
    A::T
    B::T
    n::Int
    k::Int
    char_vec::Vector{zzModRingElem}
    cache::Dict{Symbol, Any}
end

mutable struct AsymmetricProductCode <: AbstractStabilizerCodeCSS
    F::CTFieldTypes
    S1::AbstractSubsystemCode
    S2::AbstractSubsystemCode
    n::Int
    k::Int
    char_vec::Vector{zzModRingElem}
    cache::Dict{Symbol, Any}
end

mutable struct SymmetricProductCode <: AbstractStabilizerCodeCSS
    F::CTFieldTypes
    vec_S::Vector{<:AbstractSubsystemCode}
    D::Int
    n::Int
    k::Int
    char_vec::Vector{zzModRingElem}
    cache::Dict{Symbol, Any}
end

mutable struct HomologicalProductCode <: AbstractStabilizerCodeCSS
    F::CTFieldTypes
    S1::AbstractStabilizerCode
    S2::AbstractStabilizerCode
    U::CTMatrixTypes
    V::CTMatrixTypes
    n::Int
    k::Int
    char_vec::Vector{zzModRingElem}
    cache::Dict{Symbol, Any}
end

#############################
# Quantum/GeneralizedToricCode.jl
#############################


struct MonomialCode <: AbstractMonomialCode
      LR::AbstractAlgebra.Generic.LaurentMPolyWrapRing{fpFieldElem, fpMPolyRing}
      F::CTFieldTypes
      n::Int
      f::CTLRPolyElem
      g::CTLRPolyElem
end

struct FiniteMonomialCode <: AbstractMonomialCode
      LR::AbstractAlgebra.Generic.LaurentMPolyWrapRing{fpFieldElem, fpMPolyRing}
      F::CTFieldTypes
      n::Int
      k::Int
      f::CTLRPolyElem
      g::CTLRPolyElem
      a1::Tuple{Int, Int}
      a2::Tuple{Int, Int}
end

struct BivariateBicycleCode <: AbstractBivariateBicycleCode
      R::CTPolyRing
      F::CTFieldTypes
      n::Int
      k::Int
      f::CTPolyRingElem
      g::CTPolyRingElem
      l::Int
      m::Int
end

struct CoprimeBivariateBicycleCode <: AbstractBivariateBicycleCode
      R::CTPolyRing
      F::CTFieldTypes
      n::Int
      k::Int
      f::CTPolyRingElem
      g::CTPolyRingElem
      l::Int
      m::Int
end

struct GeneralizedToricCode <: AbstractGeneralizedToricCode
      LR::AbstractAlgebra.Generic.LaurentMPolyWrapRing{fpFieldElem, fpMPolyRing}
      F::CTFieldTypes
      n::Int
      f::CTLRPolyElem
      g::CTLRPolyElem
end

struct FiniteGeneralizedToricCode <: AbstractGeneralizedToricCode
      LR::AbstractAlgebra.Generic.LaurentMPolyWrapRing{fpFieldElem, fpMPolyRing}
      F::CTFieldTypes
      n::Int
      k::Int
      f::CTLRPolyElem
      g::CTLRPolyElem
      a1::Tuple{Int, Int}
      a2::Tuple{Int, Int}
end

struct Generalized3DToricCode <: AbstractGeneralized3DToricCode
      LR::AbstractAlgebra.Generic.LaurentMPolyWrapRing{fpFieldElem, fpMPolyRing}
      F::CTFieldTypes
      f::CTLRPolyElem
      g::CTLRPolyElem
end

struct FiniteGeneralized3DToricCode <: AbstractGeneralized3DToricCode
      LR::AbstractAlgebra.Generic.LaurentMPolyWrapRing{fpFieldElem, fpMPolyRing}
      F::CTFieldTypes
      f::CTLRPolyElem
      g::CTLRPolyElem
      a1::Tuple{Int, Int}
      a2::Tuple{Int, Int}
      l::Int
end

#############################
          # traits
#############################

const CSSTypes = Union{AbstractSubsystemCodeCSS, AbstractStabilizerCodeCSS, AbstractGraphStateStabilizerCSS, AbstractGraphStateSubsystemCSS, AbstractHypergraphProductCode}
const GraphStateTypes = Union{AbstractGraphStateSubsystem, AbstractGraphStateSubsystemCSS, AbstractGraphStateStabilizer, AbstractGraphStateStabilizerCSS}

"""
$(TYPEDEF)

Holy-trait function and root trait type that map a quantum code type to `HasLogicals` or `HasNoLogicals`. Dispatch on this trait instead of testing concrete code types.
"""
abstract type LogicalTrait end

"""
$(TYPEDEF)

Singleton trait returned by `LogicalTrait` for code types that carry logical operators.
"""
struct HasLogicals <: LogicalTrait end

"""
$(TYPEDEF)

Singleton trait returned by `LogicalTrait` for graph-state code types, which do not carry logical operators.
"""
struct HasNoLogicals <: LogicalTrait end
LogicalTrait(::Type{T}) where {T <: AbstractSubsystemCode} = HasLogicals()
LogicalTrait(::Type{T}) where {T <: GraphStateTypes} = HasNoLogicals()

"""
$(TYPEDEF)

Holy-trait function and root trait type that map a quantum code type to `HasGauges` or `HasNoGauges`. Dispatch on this trait instead of testing concrete code types.
"""
abstract type GaugeTrait end

"""
$(TYPEDEF)

Singleton trait returned by `GaugeTrait` for subsystem-code types with gauge operators.
"""
struct HasGauges <: GaugeTrait end

"""
$(TYPEDEF)

Singleton trait returned by `GaugeTrait` for stabilizer-code types without gauge operators.
"""
struct HasNoGauges <: GaugeTrait end
GaugeTrait(::Type{T}) where {T <: AbstractSubsystemCode} = HasGauges()
GaugeTrait(::Type{T}) where {T <: AbstractStabilizerCode} = HasNoGauges()

"""
$(TYPEDEF)

Holy-trait function and root trait type that map a quantum code type to `IsCSS` or `IsNotCSS`. Dispatch on this trait instead of testing concrete code types.
"""
abstract type CSSTrait end

"""
$(TYPEDEF)

Singleton trait returned by `CSSTrait` for code types in `CSSTypes`.
"""
struct IsCSS <: CSSTrait end

"""
$(TYPEDEF)

Singleton trait returned by `CSSTrait` for non-CSS subsystem-code types.
"""
struct IsNotCSS <: CSSTrait end
CSSTrait(::Type{T}) where {T <: AbstractSubsystemCode} = IsNotCSS()
CSSTrait(::Type{T}) where {T <: CSSTypes} = IsCSS()
