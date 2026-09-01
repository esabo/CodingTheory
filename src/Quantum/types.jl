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

abstract type AbstractGeneralized3DToricCode <: AbstractStabilizerCodeCSS end

# AbstractQuantumLDPCCode, AbstractQuantumLDPCCSSCode?

abstract type AbstractQuantumNoiseChannel <: AbstractChannel end

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
# Quantum/product_codes.jl
#############################


#############################
# Quantum/GeneralizedToricCode.jl
#############################

mutable struct BivariateBicycleCode{T} <: AbstractStabilizerCodeCSS
    LR::T # LaurentMPolyWrapRing
    F::CTFieldTypes
    a::Any # CTLRPolyElem
    b::Any # CTLRPolyElem
    a1::Tuple{Int, Int}
    a2::Tuple{Int, Int}
    n::Int
    k::Union{Int, Missing}
    d::Union{Int, Missing}
    l_bound::Int
    u_bound::Int
    cache::Dict{Symbol, Any}
end

mutable struct Generalized3DToricCode{T} <: AbstractGeneralized3DToricCode
    LR::T # LaurentMPolyWrapRing
    F::CTFieldTypes
    a::Any # CTLRPolyElem
    b::Any # CTLRPolyElem
    a1::Tuple{Int, Int}
    a2::Tuple{Int, Int}
    l_z::Int
    n::Int
    k::Union{Int, Missing}
    d::Union{Int, Missing}
    l_bound::Int
    u_bound::Int
    cache::Dict{Symbol, Any}
end

mutable struct QuantumConcatenatedCode <: AbstractStabilizerCode
    outer_code::AbstractStabilizerCode
    inner_code::AbstractStabilizerCode
    n::Int
    k::Int
    d::Union{Int, Missing}
    l_bound::Int
    u_bound::Int
    cache::Dict{Symbol, Any}
end

mutable struct GaugeFixedCode <: AbstractStabilizerCode
    subsystem_code::AbstractSubsystemCode
    choice::Symbol
    n::Int
    k::Int
    d::Union{Int, Missing}
    l_bound::Int
    u_bound::Int
    cache::Dict{Symbol, Any}
end

mutable struct HypergraphProductCode <: AbstractStabilizerCodeCSS
    C1::AbstractLinearCode
    C2::AbstractLinearCode
    n::Int
    k::Int
    d::Union{Int, Missing}
    l_bound::Int
    u_bound::Int
    cache::Dict{Symbol, Any}
end

mutable struct GeneralizedShorCode <: AbstractSubsystemCode
    C1::AbstractLinearCode
    C2::AbstractLinearCode
    n::Int
    k::Int
    r::Union{Int, Missing} # Gauge qubits, computed lazily
    d::Union{Int, Missing}
    l_bound::Int
    u_bound::Int
    cache::Dict{Symbol, Any}
end

mutable struct HyperBicycleCodeCSS{T <: CTMatrixTypes} <: AbstractStabilizerCodeCSS
    a::Vector{T}
    b::Vector{T}
    χ::Int
    n::Int
    k::Union{Int, Missing}
    d::Union{Int, Missing}
    l_bound::Int
    u_bound::Int
    cache::Dict{Symbol, Any}
end

mutable struct HyperBicycleCode{T <: CTMatrixTypes} <: AbstractStabilizerCode
    a::Vector{T}
    b::Vector{T}
    χ::Int
    n::Int
    k::Union{Int, Missing}
    d::Union{Int, Missing}
    l_bound::Int
    u_bound::Int
    cache::Dict{Symbol, Any}
end

# Intercept property access to compute 'k' lazily when requested
for T in (:HyperBicycleCodeCSS, :HyperBicycleCode)
    @eval begin
        function Base.getproperty(S::$T, prop::Symbol)
            if prop == :k
                k_val = getfield(S, :k)
                if ismissing(k_val)
                    # For CSS, k = n - rank(H_X) - rank(H_Z)
                    # For non-CSS Stabilizer, k = n - rank(stabs)
                    if $T == HyperBicycleCodeCSS
                        H_X = X_stabilizers(S)
                        H_Z = Z_stabilizers(S)
                        k_val = S.n - rank(H_X) - rank(H_Z)
                    else
                        stabs = stabilizers(S)
                        k_val = S.n - rank(stabs)
                    end
                    setfield!(S, :k, k_val)
                    return k_val
                end
                return k_val
            else
                return getfield(S, prop)
            end
        end
    end
end

mutable struct GeneralizedBicycleCode{T <: CTMatrixTypes} <: AbstractStabilizerCodeCSS
    A::T
    B::T
    n::Int
    k::Union{Int, Missing}
    d::Union{Int, Missing}
    l_bound::Int
    u_bound::Int
    cache::Dict{Symbol, Any}
end

# Intercept property access to compute 'k' lazily when requested
function Base.getproperty(S::GeneralizedBicycleCode, prop::Symbol)
    if prop == :k
        k_val = getfield(S, :k)
        if ismissing(k_val)
            H_X = X_stabilizers(S)
            H_Z = Z_stabilizers(S)
            k_val = S.n - rank(H_X) - rank(H_Z)
            setfield!(S, :k, k_val)
            return k_val
        end
        return k_val
    else
        return getfield(S, prop)
    end
end

mutable struct LiftedProductCode{T} <: AbstractStabilizerCodeCSS
    A::T
    B::T
    n::Int
    k::Union{Int, Missing}
    d::Union{Int, Missing}
    l_bound::Int
    u_bound::Int
    cache::Dict{Symbol, Any}
end

mutable struct BiasTailoredLiftedProductCode{T} <: AbstractStabilizerCode
    A::T
    B::T
    n::Int
    k::Union{Int, Missing}
    d::Union{Int, Missing}
    l_bound::Int
    u_bound::Int
    cache::Dict{Symbol, Any}
end

for T in (:LiftedProductCode, :BiasTailoredLiftedProductCode)
    @eval begin
        function Base.getproperty(S::$T, prop::Symbol)
            if prop == :k
                k_val = getfield(S, :k)
                if ismissing(k_val)
                    if $T == LiftedProductCode
                        H_X = X_stabilizers(S)
                        H_Z = Z_stabilizers(S)
                        k_val = S.n - rank(H_X) - rank(H_Z)
                    else
                        stabs = stabilizers(S)
                        k_val = S.n - rank(stabs)
                    end
                    setfield!(S, :k, k_val)
                    return k_val
                end
                return k_val
            else
                return getfield(S, prop)
            end
        end
    end
end

mutable struct AsymmetricProductCode <: AbstractStabilizerCodeCSS
    S1::AbstractSubsystemCode
    S2::AbstractSubsystemCode
    n::Int
    k::Union{Int, Missing}
    d::Union{Int, Missing}
    l_bound::Int
    u_bound::Int
    cache::Dict{Symbol, Any}
end

mutable struct SymmetricProductCode <: AbstractStabilizerCodeCSS
    vec_S::Vector{<:AbstractSubsystemCode}
    D::Int
    n::Int
    k::Union{Int, Missing}
    d::Union{Int, Missing}
    l_bound::Int
    u_bound::Int
    cache::Dict{Symbol, Any}
end

for T in (:AsymmetricProductCode, :SymmetricProductCode)
    @eval begin
        function Base.getproperty(S::$T, prop::Symbol)
            if prop == :k
                k_val = getfield(S, :k)
                if ismissing(k_val)
                    H_X = X_stabilizers(S)
                    H_Z = Z_stabilizers(S)
                    k_val = S.n - rank(H_X) - rank(H_Z)
                    setfield!(S, :k, k_val)
                    return k_val
                end
                return k_val
            else
                return getfield(S, prop)
            end
        end
    end
end

mutable struct HomologicalProductCode <: AbstractStabilizerCodeCSS
    S1::AbstractStabilizerCode
    S2::AbstractStabilizerCode
    U::CTMatrixTypes
    V::CTMatrixTypes
    n::Int
    k::Union{Int, Missing}
    d::Union{Int, Missing}
    l_bound::Int
    u_bound::Int
    cache::Dict{Symbol, Any}
end

# Intercept property access to compute 'k' lazily when requested
function Base.getproperty(S::HomologicalProductCode, prop::Symbol)
    if prop == :k
        k_val = getfield(S, :k)
        if ismissing(k_val)
            H_X = X_stabilizers(S)
            H_Z = Z_stabilizers(S)
            k_val = S.n - rank(H_X) - rank(H_Z)
            setfield!(S, :k, k_val)
            return k_val
        end
        return k_val
    else
        return getfield(S, prop)
    end
end

mutable struct CoprimeBivariateBicycleCode <: AbstractStabilizerCodeCSS
    R::CTPolyRing
    F::CTFieldTypes
    a::CTPolyRingElem
    b::CTPolyRingElem
    N::Int
    n::Int
    k::Union{Int, Missing}
    d::Union{Int, Missing}
    l_bound::Int
    u_bound::Int
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
