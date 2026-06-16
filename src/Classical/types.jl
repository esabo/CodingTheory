# Copyright (c) 2023 - 2026 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
      # abstract types
#############################

abstract type AbstractCode end
abstract type AbstractNonadditiveCode <: AbstractCode end
abstract type AbstractNonlinearCode <: AbstractCode end
abstract type AbstractAdditiveCode <: AbstractCode end
abstract type AbstractLinearCode <: AbstractAdditiveCode end
abstract type AbstractMatrixProductCode <: AbstractLinearCode end
abstract type AbstractReedMullerCode <: AbstractLinearCode end
abstract type AbstractCyclicCode <: AbstractLinearCode end
abstract type AbstractBCHCode <: AbstractCyclicCode end
abstract type AbstractReedSolomonCode <: AbstractBCHCode end
abstract type AbstractCyclicCode2D <: AbstractCyclicCode end
abstract type AbstractQuasiCyclicCode <: AbstractLinearCode end
abstract type AbstractGeneralizedReedSolomonCode <: AbstractLinearCode end
abstract type AbstractAlgebraicGeometryCode <: AbstractLinearCode end
abstract type AbstractConcatenatedCode <: AbstractLinearCode end
abstract type AbstractAlternateCode <: AbstractLinearCode end
abstract type AbstractGoppaCode <: AbstractAlternateCode end
abstract type AbstractGeneralizedSrivastavaCode <: AbstractAlternateCode end
abstract type AbstractTwistedReedSolomonCode <: AbstractLinearCode end
abstract type AbstractConcatenatedCode <: AbstractLinearCode end
abstract type AbstractTannerCode <: AbstractLinearCode end

#############################
      # concrete types
#############################

function Base.getproperty(C::AbstractLinearCode, sym::Symbol)
    # 1. Physical fields
    if sym in fieldnames(typeof(C))
        return getfield(C, sym)
    end
    
    # 2. Check Cache
    cache = getfield(C, :cache)
    if haskey(cache, sym)
        return cache[sym]
    end
    
    # 3. Lazy Evaluation Routing
    if sym == :G
        return generator_matrix(C)
    elseif sym == :H
        return parity_check_matrix(C)
    elseif sym == :weight_enum
        return weight_enumerator(C)
    elseif sym == :weight_dist
        return weight_distribution(C)
    end
    
    throw(ErrorException("type $(typeof(C)) has no field $sym and it is not in the cache."))
end

# We also overload setproperty! so C.G_stand = ... routes to the cache
function Base.setproperty!(C::AbstractLinearCode, sym::Symbol, val)
    if sym in fieldnames(typeof(C))
        setfield!(C, sym, val)
    else
        getfield(C, :cache)[sym] = val
    end
end

#############################
       # linearcode.jl
#############################

# struct WeightEnumerator
#       polynomial::Union{ZZMPolyRingElem, Nemo.AbsSimpleNumFieldElem}
#       type::Symbol
# end
    
struct HammingWeightEnumerator
    n::Int
    counts::Dict{Int, BigInt}
end

struct CompleteWeightEnumerator
    n::Int
    q::Int
    elements::Vector{Any}
    counts::Dict{Tuple, BigInt}
end

# mutable struct ExtendedQRCode <: AbstractLinearCode
#     F::CTFieldTypes # base field
#     n::Int # length
#     k::Int # dimension
#     d::Union{Int, Missing} # minimum distance
#     l_bound::Int # lower bound on d
#     u_bound::Int # upper bound on d
#     G::CTMatrixTypes
#     H::CTMatrixTypes
#     G_stand::CTMatrixTypes
#     H_stand::CTMatrixTypes
#     P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
#     weight_enum::Union{WeightEnumerator, Missing}
# end

mutable struct ExtendedQRCode <: AbstractLinearCode
    F::CTFieldTypes # base field
    n::Int # length
    k::Int # dimension
    d::Union{Int, Missing} # minimum distance
    l_bound::Int # lower bound on d
    u_bound::Int # upper bound on d
    G::CTMatrixTypes
    H::CTMatrixTypes
    cache::Dict{Symbol, Any}
end

# mutable struct ProductCode <: AbstractLinearCode
#     F::CTFieldTypes # base field
#     n::Int # length
#     k::Int # dimension
#     d::Union{Int, Missing} # minimum distance
#     l_bound::Int # lower bound on d
#     u_bound::Int # upper bound on d
#     G::CTMatrixTypes
#     H::CTMatrixTypes
#     G_stand::CTMatrixTypes
#     H_stand::CTMatrixTypes
#     P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
#     weight_enum::Union{WeightEnumerator, Missing}
# end

mutable struct ProductCode <: AbstractLinearCode
    C1::AbstractLinearCode # The row code
    C2::AbstractLinearCode # The column code
    F::CTFieldTypes 
    n::Int 
    k::Int 
    d::Union{Int, Missing}
    l_bound::Int 
    u_bound::Int 
    cache::Dict{Symbol, Any}
end

# mutable struct LinearCode <: AbstractLinearCode
#     F::CTFieldTypes # base field
#     n::Int # length
#     k::Int # dimension
#     d::Union{Int, Missing} # minimum distance
#     l_bound::Int # lower bound on d
#     u_bound::Int # upper bound on d
#     G::CTMatrixTypes
#     H::CTMatrixTypes
#     G_stand::CTMatrixTypes
#     H_stand::CTMatrixTypes
#     P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
#     weight_enum::Union{WeightEnumerator, Missing}
# end

mutable struct LinearCode <: AbstractLinearCode
    F::CTFieldTypes 
    n::Int 
    k::Int 
    d::Union{Int, Missing} 
    l_bound::Int 
    u_bound::Int 
    cache::Dict{Symbol, Any}
end

#############################
    # MatrixProductCode.jl
#############################

# mutable struct MatrixProductCode <: AbstractMatrixProductCode
#     F::CTFieldTypes # base field
#     n::Int # length
#     k::Int # dimension
#     d::Union{Int, Missing} # minimum distance
#     l_bound::Int # lower bound on d
#     u_bound::Int # upper bound on d
#     G::CTMatrixTypes
#     H::CTMatrixTypes
#     G_stand::CTMatrixTypes
#     H_stand::CTMatrixTypes
#     P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
#     weight_enum::Union{WeightEnumerator, Missing}
#     Cvec::Vector{AbstractLinearCode}
#     A::fqPolyRepMatrix
# end

#############################
       # ReedMuller.jl
#############################

# mutable struct ReedMullerCode <: AbstractReedMullerCode
#     F::CTFieldTypes
#     n::Int # length
#     k::Int # dimension
#     d::Union{Int, Missing} # minimum distance
#     l_bound::Int # lower bound on d
#     u_bound::Int # upper bound on d
#     r::Integer # order
#     m::Integer # number of variables
#     G::CTMatrixTypes
#     H::CTMatrixTypes
#     G_stand::CTMatrixTypes
#     H_stand::CTMatrixTypes
#     P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
#     weight_enum::Union{WeightEnumerator, Missing}
# end
  
#############################
       # cycliccode.jl
#############################

# mutable struct CyclicCode <: AbstractCyclicCode
#     F::CTFieldTypes # base field
#     E::CTFieldTypes # splitting field
#     R::CTPolyRing # polynomial ring of generator polynomial
#     β::CTFieldElem # n-th root of primitive element of splitting field
#     n::Int # length
#     k::Int # dimension
#     d::Union{Int, Missing} # minimum distance
#     b::Int # offset
#     δ::Int # BCH bound
#     HT::Int # Hartmann-Tzeng refinement
#     l_bound::Int # lower bound on d
#     u_bound::Int # upper bound on d
#     qcosets::Vector{Vector{Int}}
#     qcosets_reps::Vector{Int}
#     def_set::Vector{Int}
#     g::CTPolyRingElem
#     h::CTPolyRingElem
#     e::CTPolyRingElem
#     G::CTMatrixTypes
#     H::CTMatrixTypes
#     G_stand::CTMatrixTypes
#     H_stand::CTMatrixTypes
#     P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
#     weight_enum::Union{WeightEnumerator, Missing}
# end

mutable struct CyclicCode <: AbstractCyclicCode
    F::CTFieldTypes
    E::CTFieldTypes
    R::CTPolyRing
    β::CTFieldElem
    n::Int
    k::Int
    l_bound::Int
    u_bound::Int
    qcosets::Vector{Vector{Int}}
    qcosets_reps::Vector{Int}
    def_set::Vector{Int}
    g::CTPolyRingElem
    h::CTPolyRingElem
    e::CTPolyRingElem
    cache::Dict{Symbol, Any}
end

# mutable struct BCHCode <: AbstractBCHCode
#     F::CTFieldTypes # base field
#     E::CTFieldTypes # splitting field
#     R::CTPolyRing # polynomial ring of generator polynomial
#     β::CTFieldElem # n-th root of primitive element of splitting field
#     n::Int # length
#     k::Int # dimension
#     d::Union{Int, Missing} # minimum distance
#     b::Int # offset
#     δ::Int # BCH bound
#     HT::Int # Hartmann-Tzeng refinement
#     l_bound::Int # lower bound on d
#     u_bound::Int # upper bound on d
#     qcosets::Vector{Vector{Int}}
#     qcosets_reps::Vector{Int}
#     def_set::Vector{Int}
#     g::CTPolyRingElem
#     h::CTPolyRingElem
#     e::CTPolyRingElem
#     G::CTMatrixTypes
#     H::CTMatrixTypes
#     G_stand::CTMatrixTypes
#     H_stand::CTMatrixTypes
#     P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
#     weight_enum::Union{WeightEnumerator, Missing}
# end

mutable struct BCHCode <: AbstractBCHCode
    F::CTFieldTypes
    E::CTFieldTypes
    R::CTPolyRing
    β::CTFieldElem
    n::Int
    k::Int
    l_bound::Int
    u_bound::Int
    qcosets::Vector{Vector{Int}}
    qcosets_reps::Vector{Int}
    def_set::Vector{Int}
    g::CTPolyRingElem
    h::CTPolyRingElem
    e::CTPolyRingElem
    cache::Dict{Symbol, Any}
end

# mutable struct ReedSolomonCode <: AbstractReedSolomonCode
#     F::CTFieldTypes # base field
#     E::CTFieldTypes # splitting field
#     R::CTPolyRing # polynomial ring of generator polynomial
#     β::CTFieldElem # n-th root of primitive element of splitting field
#     n::Int # length
#     k::Int # dimension
#     d::Union{Int, Missing} # minimum distance
#     b::Int # offset
#     δ::Int # BCH bound
#     HT::Int # Hartmann-Tzeng refinement
#     l_bound::Int # lower bound on d
#     u_bound::Int # upper bound on d
#     qcosets::Vector{Vector{Int}}
#     qcosets_reps::Vector{Int}
#     def_set::Vector{Int}
#     g::CTPolyRingElem
#     h::CTPolyRingElem
#     e::CTPolyRingElem
#     G::CTMatrixTypes
#     H::CTMatrixTypes
#     G_stand::CTMatrixTypes
#     H_stand::CTMatrixTypes
#     P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
#     weight_enum::Union{WeightEnumerator, Missing}
# end

mutable struct ReedSolomonCode <: AbstractReedSolomonCode
    F::CTFieldTypes
    E::CTFieldTypes
    R::CTPolyRing
    β::CTFieldElem
    n::Int
    k::Int
    l_bound::Int
    u_bound::Int
    qcosets::Vector{Vector{Int}}
    qcosets_reps::Vector{Int}
    def_set::Vector{Int}
    g::CTPolyRingElem
    h::CTPolyRingElem
    e::CTPolyRingElem
    cache::Dict{Symbol, Any}
end

#############################
    # quasicycliccode.jl
#############################

# mutable struct QuasiCyclicCode <: AbstractQuasiCyclicCode
#     F::CTFieldTypes # base field
#     R::EuclideanRingResidueRing
#     n::Int # length
#     k::Int # dimension
#     d::Union{Int, Missing} # minimum distance
#     l_bound::Int # lower bound on d
#     u_bound::Int # upper bound on d
#     G::Union{CTMatrixTypes, Missing}
#     H::Union{CTMatrixTypes, Missing}
#     G_stand::Union{CTMatrixTypes, Missing}
#     H_stand::Union{CTMatrixTypes, Missing}
#     P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
#     weight_enum::Union{WeightEnumerator, Missing}
#     l::Int
#     m::Int
#     A::MatElem{<:ResElem}
#     A_type::Symbol
#     W::Matrix{Int}
#     type::Int
# end

mutable struct QuasiCyclicCode <: AbstractQuasiCyclicCode
    F::CTFieldTypes
    R::CTPolyRing
    n::Int
    k::Int
    d::Union{Int, Missing}
    l_bound::Int
    u_bound::Int
    l::Int
    m::Int
    A::CTMatrixTypes
    A_type::Symbol
    cache::Dict{Symbol, Any}
end

#############################
      # GRS_alternate.jl
#############################

# mutable struct GeneralizedReedSolomonCode <: AbstractGeneralizedReedSolomonCode
#     F::CTFieldTypes # base field
#     n::Int # length
#     k::Int # dimension
#     d::Union{Int, Missing} # minimum distance
#     l_bound::Int # lower bound on d
#     u_bound::Int # upper bound on d
#     scalars::Vector{<:CTFieldElem}
#     dual_scalars::Vector{<:CTFieldElem}
#     eval_pts::Vector{<:CTFieldElem}
#     G::CTMatrixTypes
#     H::CTMatrixTypes
#     G_stand::CTMatrixTypes
#     H_stand::CTMatrixTypes
#     P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
#     weight_enum::Union{WeightEnumerator, Missing} # TODO: should never be missing? is completely known for MDS?
# end

mutable struct GeneralizedReedSolomonCode <: AbstractGeneralizedReedSolomonCode
    F::CTFieldTypes 
    n::Int 
    k::Int 
    d::Union{Int, Missing} 
    l_bound::Int 
    u_bound::Int 
    scalars::Vector{<:CTFieldElem}
    dual_scalars::Vector{<:CTFieldElem}
    eval_pts::Vector{<:CTFieldElem}
    cache::Dict{Symbol, Any}
end

# mutable struct AlternateCode <: AbstractAlternateCode
#     F::CTFieldTypes # base field
#     E::CTFieldTypes # extension field
#     n::Int # length
#     k::Int # dimension
#     d::Union{Int, Missing} # minimum distance
#     l_bound::Int # lower bound on d
#     u_bound::Int # upper bound on d
#     scalars::Vector{<:CTFieldElem}
#     eval_pts::Vector{<:CTFieldElem}
#     G::CTMatrixTypes
#     H::CTMatrixTypes
#     G_stand::CTMatrixTypes
#     H_stand::CTMatrixTypes
#     P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
#     weight_enum::Union{WeightEnumerator, Missing}
# end

mutable struct AlternateCode <: AbstractAlternateCode
    F::CTFieldTypes 
    E::CTFieldTypes 
    n::Int 
    k::Int 
    d::Union{Int, Missing} 
    l_bound::Int 
    u_bound::Int 
    scalars::Vector{<:CTFieldElem}
    eval_pts::Vector{<:CTFieldElem}
    cache::Dict{Symbol, Any}
end

# mutable struct GeneralizedSrivastavaCode <: AbstractGeneralizedSrivastavaCode
#     F::CTFieldTypes # base field
#     E::CTFieldTypes # extension field
#     n::Int # length
#     k::Int # dimension
#     d::Union{Int, Missing} # minimum distance
#     l_bound::Int # lower bound on d
#     u_bound::Int # upper bound on d
#     a::Vector{<:CTFieldElem}
#     w::Vector{<:CTFieldElem}
#     z::Vector{<:CTFieldElem}
#     t::Int
#     G::CTMatrixTypes
#     H::CTMatrixTypes
#     G_stand::CTMatrixTypes
#     H_stand::CTMatrixTypes
#     P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
#     weight_enum::Union{WeightEnumerator, Missing}
# end

mutable struct GeneralizedSrivastavaCode <: AbstractGeneralizedSrivastavaCode
    F::CTFieldTypes 
    E::CTFieldTypes 
    n::Int 
    k::Int 
    d::Union{Int, Missing} 
    l_bound::Int 
    u_bound::Int 
    a::Vector{<:CTFieldElem}
    w::Vector{<:CTFieldElem}
    z::Vector{<:CTFieldElem}
    t::Int
    cache::Dict{Symbol, Any}
end

#############################
         # Goppa.jl
#############################

# mutable struct GoppaCode <: AbstractGoppaCode
#     F::CTFieldTypes # base field
#     E::CTFieldTypes # extension field
#     n::Int # length
#     k::Int # dimension
#     d::Union{Int, Missing} # minimum distance
#     l_bound::Int # lower bound on d
#     u_bound::Int # upper bound on d
#     G::CTMatrixTypes
#     H::CTMatrixTypes
#     G_stand::CTMatrixTypes
#     H_stand::CTMatrixTypes
#     P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
#     weight_enum::Union{WeightEnumerator, Missing}
#     L::Vector{CTFieldElem}
#     g::FqPolyRingElem
# end

#############################
     # concatenation.jl
#############################

# mutable struct ConcatenatedCode <: AbstractLinearCode
#     C_in::Union{AbstractLinearCode, Vector{<:AbstractLinearCode}}
#     C_out::Union{AbstractLinearCode, Vector{<:AbstractLinearCode}}
#     type::Union{Symbol, Vector{Symbol}}
#     basis::Union{Missing, Vector{Union{Missing, <:CTFieldElem, Vector{<:CTFieldElem}}}}
#     dual_basis::Union{Missing, Vector{Union{Missing, <:CTFieldElem, Vector{<:CTFieldElem}}}}
#     F::CTFieldTypes # base field
#     n::Int # length
#     k::Int # dimension
#     d::Union{Int, Missing} # minimum distance
#     l_bound::Int # lower bound on d
#     u_bound::Int # upper bound on d
#     G::CTMatrixTypes
#     H::CTMatrixTypes
#     G_stand::CTMatrixTypes
#     H_stand::CTMatrixTypes
#     P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
#     weight_enum::Union{WeightEnumerator, Missing}
# end

mutable struct ConcatenatedCode <: AbstractLinearCode
    C_in::Union{AbstractLinearCode, Vector{<:AbstractLinearCode}}
    C_out::Union{AbstractLinearCode, Vector{<:AbstractLinearCode}}
    type::Union{Symbol, Vector{Symbol}}
    basis::Union{Missing, Vector{Union{Missing, <:CTFieldElem, Vector{<:CTFieldElem}}}}
    dual_basis::Union{Missing, Vector{Union{Missing, <:CTFieldElem, Vector{<:CTFieldElem}}}}
    F::CTFieldTypes 
    n::Int 
    k::Int 
    d::Union{Int, Missing} 
    l_bound::Int 
    u_bound::Int 
    cache::Dict{Symbol, Any}
end

#############################
   # TwistedReedSolomon.jl
#############################

# mutable struct TwistedReedSolomonCode <: AbstractTwistedReedSolomonCode
#     F::CTFieldTypes # base field
#     n::Int # length
#     k::Int # dimension
#     d::Union{Int, Missing} # minimum distance
#     l_bound::Int # lower bound on d
#     u_bound::Int # upper bound on d
#     G::CTMatrixTypes
#     H::CTMatrixTypes
#     G_stand::CTMatrixTypes
#     H_stand::CTMatrixTypes
#     P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
#     weight_enum::Union{WeightEnumerator, Missing}
#     α::Vector{T} where T <: CTFieldElem
#     t::Vector{Int}
#     h::Vector{Int}
#     η::Vector{T} where T <: CTFieldElem
#     l::Int
# end

mutable struct TwistedReedSolomonCode <: AbstractTwistedReedSolomonCode
    F::CTFieldTypes 
    n::Int 
    k::Int 
    d::Union{Int, Missing} 
    l_bound::Int 
    u_bound::Int
    α::Vector{T} where T <: CTFieldElem
    t::Vector{Int}
    h::Vector{Int}
    η::Vector{T} where T <: CTFieldElem
    l::Int
    cache::Dict{Symbol, Any}
end

# ==============================================================================
# GEOMETRIC & SPORADIC CODE TYPES
# ==============================================================================

mutable struct HammingCode <: AbstractLinearCode
    F::CTFieldTypes 
    n::Int 
    k::Int 
    d::Union{Int, Missing} 
    l_bound::Int 
    u_bound::Int 
    r::Int             # Preserved structural parameter
    cache::Dict{Symbol, Any}
end

mutable struct SimplexCode <: AbstractLinearCode
    F::CTFieldTypes 
    n::Int 
    k::Int 
    d::Union{Int, Missing} 
    l_bound::Int 
    u_bound::Int 
    r::Int             # Preserved structural parameter
    cache::Dict{Symbol, Any}
end

mutable struct MacDonaldCode <: AbstractLinearCode
    F::CTFieldTypes 
    n::Int 
    k::Int 
    d::Union{Int, Missing} 
    l_bound::Int 
    u_bound::Int 
    u::Int             # The dimension of the punctured subcode
    cache::Dict{Symbol, Any}
end

# ==============================================================================
# COMPOSITE CODE TYPES
# ==============================================================================

mutable struct PlotkinCode <: AbstractLinearCode
    C1::AbstractLinearCode
    C2::AbstractLinearCode
    F::CTFieldTypes 
    n::Int 
    k::Int 
    d::Union{Int, Missing} 
    l_bound::Int 
    u_bound::Int 
    cache::Dict{Symbol, Any}
end

mutable struct DirectSumCode <: AbstractLinearCode
    C1::AbstractLinearCode
    C2::AbstractLinearCode
    F::CTFieldTypes 
    n::Int 
    k::Int 
    d::Union{Int, Missing} 
    l_bound::Int 
    u_bound::Int 
    cache::Dict{Symbol, Any}
end

mutable struct TensorProductCode <: AbstractLinearCode
    C1::AbstractLinearCode
    C2::AbstractLinearCode
    F::CTFieldTypes 
    n::Int 
    k::Int 
    d::Union{Int, Missing} 
    l_bound::Int 
    u_bound::Int 
    cache::Dict{Symbol, Any}
end

mutable struct MultilevelConcatenatedCode <: AbstractConcatenatedCode
    C_outs::Vector{<:AbstractLinearCode}
    C_ins::Vector{<:AbstractLinearCode}
    types::Vector{Symbol}
    bases::Vector{Union{Vector{<:CTFieldElem}, Missing}}
    dual_bases::Vector{Union{Vector{<:CTFieldElem}, Missing}}
    F::CTFieldTypes
    n::Int
    k::Int
    d::Union{Int, Missing}
    l_bound::Int
    u_bound::Int
    cache::Dict{Symbol, Any}
end

mutable struct GabidulinCode <: AbstractLinearCode
    F::CTFieldTypes
    E::CTFieldTypes
    n::Int
    k::Int
    d::Union{Int, Missing}
    l_bound::Int
    u_bound::Int
    eval_pts::Vector{<:CTFieldElem}
    s::Int
    cache::Dict{Symbol, Any}
end

mutable struct GoppaCode <: AbstractGoppaCode
    F::CTFieldTypes
    E::CTFieldTypes
    n::Int
    k::Int
    d::Union{Int, Missing}
    l_bound::Int
    u_bound::Int
    L::Vector{<:CTFieldElem}
    g::CTPolyRingElem
    cache::Dict{Symbol, Any}
end

mutable struct MatrixProductCode <: AbstractLinearCode
    C::Vector{<:AbstractLinearCode}
    A::CTMatrixTypes
    F::CTFieldTypes
    n::Int
    k::Int
    d::Union{Int, Missing}
    l_bound::Int
    u_bound::Int
    cache::Dict{Symbol, Any}
end

mutable struct ReedMullerCode <: AbstractLinearCode
    F::CTFieldTypes
    n::Int
    k::Int
    d::Union{Int, Missing}
    l_bound::Int
    u_bound::Int
    r::Int
    m::Int
    cache::Dict{Symbol, Any}
end

mutable struct TannerCode <: AbstractTannerCode
    F::CTFieldTypes
    n::Int
    k::Int
    d::Union{Int, Missing}
    l_bound::Int
    u_bound::Int
    cache::Dict{Symbol, Any}
end
