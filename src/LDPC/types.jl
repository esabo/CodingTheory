# Copyright (c) 2023 - 2024 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
      # abstract types
#############################

abstract type AbstractLDPCCode <: AbstractLinearCode end
# abstract type AbstractNoiseChannel end
# abstract type AbstractClassicalNoiseChannel <: AbstractNoiseChannel end
# abstract type AbstractBinaryErasureChannel <: AbstractClassicalNoiseChannel end
# abstract type AbstractBinarySymmetricChannel <: AbstractClassicalNoiseChannel end
# abstract type AbstractBAWGNChannel <: AbstractClassicalNoiseChannel end
abstract type AbstractChannel end
abstract type AbstractDiscreteChannel <: AbstractChannel end
abstract type AbstractContinuousChannel <: AbstractChannel end

#############################
      # concrete types
#############################

#############################
        # LDPC/codes.jl
#############################

mutable struct LDPCCode <: AbstractLDPCCode
    F::CTFieldTypes
    n::Int
    k::Union{Int, Missing} # Exact dimension (lazy)
    k_design::Int          # Design dimension (n - m)
    d::Union{Int, Missing}
    l_bound::Int
    u_bound::Int
    H::CTMatrixTypes 
    λ::QQPolyRingElem
    ρ::QQPolyRingElem
    is_reg::Bool
    cache::Dict{Symbol, Any}
end

#############################
      # LDPC/channels.jl
#############################

"""
$(TYPEDSIGNATURES)

Binary Erasure Channel (BEC) with erasure probability `ε`.
"""
struct BinaryErasureChannel <: AbstractDiscreteChannel
    ε::Float64
    function BinaryErasureChannel(ε::Float64)
        0.0 <= ε <= 1.0 || throw(DomainError("Erasure probability must be in [0, 1]"))
        new(ε)
    end
end
const BEC = BinaryErasureChannel

"""
$(TYPEDSIGNATURES)

Binary Symmetric Channel (BSC) with crossover probability `p`.
"""
struct BinarySymmetricChannel <: AbstractDiscreteChannel
    p::Float64
    function BinarySymmetricChannel(p::Float64)
        0.0 <= p <= 1.0 || throw(DomainError("Crossover probability must be in [0, 1]"))
        new(p)
    end
end
const BSC = BinarySymmetricChannel

"""
$(TYPEDSIGNATURES)

Binary-Input Additive White Gaussian Noise (BIAWGN) Channel.
Defined by noise standard deviation `σ`.
"""
struct BAWGNChannel <: AbstractContinuousChannel
    σ::Float64
    function BAWGNChannel(σ::Float64)
        σ > 0.0 || throw(DomainError("Standard deviation must be strictly positive"))
        new(σ)
    end
end
const BAWGNC = BAWGNChannel

struct ZChannel <: AbstractDiscreteChannel
    param::Float64
    function ZChannel(p::Float64)
        0.0 <= p <= 1.0 || throw(DomainError(p, "Transition probability p must be in [0, 1]"))
        new(p)
    end
end

struct RayleighFadingChannel <: AbstractContinuousChannel
    param::Float64 # Base noise standard deviation σ (assuming E[a^2] = 1)
    function RayleighFadingChannel(sigma::Float64)
        sigma > 0.0 || throw(DomainError(sigma, "Noise standard deviation must be positive"))
        new(sigma)
    end
end

#############################
    # LDPC/ensembles.jl
#############################

mutable struct LDPCEnsemble
    λ::PolyRingElem
    ρ::PolyRingElem
    L::PolyRingElem
    R::PolyRingElem
    l_avg::Float64
    r_avg::Float64
    design_rate::Float64
    density_evo::Dict{AbstractChannel, NTuple{2, Vector{Float64}}}
    threshold::Dict{Type, Float64}
end

"""
A structural representation of a Multi-Edge Type (MET) LDPC Ensemble.
"""
struct METEnsemble
    num_edge_types::Int
    
    # Each row is a Variable Node class.
    # Col 1: Fraction of total nodes belonging to this class
    # Col 2: Channel LLR multiplier (1.0 = transmitted, 0.0 = punctured)
    # Col 3 to end: Number of connections to each edge type
    var_profiles::Matrix{Float64}
    
    # Each row is a Check Node class.
    # Col 1: Fraction of total nodes belonging to this class
    # Col 2 to end: Number of connections to each edge type
    chk_profiles::Matrix{Float64}
end

# Define a Union type so we don't duplicate dispatches for LDPC and MET ensembles
const AbstractLDPCFamily = Union{LDPCEnsemble, METEnsemble}
