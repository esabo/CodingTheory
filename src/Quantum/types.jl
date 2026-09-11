# Copyright (c) 2023 - 2024 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

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
An algebraic three-dimensional generalized toric-code datum.
"""
struct Generalized3DToricCode{T, U, V}
    LR::T
    F::CTFieldTypes
    a::U
    b::V
end

"""
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

mutable struct LiftedProductCode{T} <: AbstractStabilizerCodeCSS
    F::CTFieldTypes
    A::T
    B::T
    n::Int
    k::Int
    char_vec::Vector{zzModRingElem}
    cache::Dict{Symbol, Any}
end

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
