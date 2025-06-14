# Copyright (c) 2023 - 2024 Eric Sabo, Benjamin Ide
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
abstract type AbstractQuasiCyclicCode <: AbstractLinearCode end
abstract type AbstractGeneralizedReedSolomonCode <: AbstractLinearCode end
abstract type AbstractAlgebraicGeometryCode <: AbstractLinearCode end
abstract type AbstractConcatenatedCode <: AbstractLinearCode end
abstract type AbstractAlternateCode <: AbstractLinearCode end
abstract type AbstractGoppaCode <: AbstractAlternateCode end
abstract type AbstractGeneralizedSrivastavaCode <: AbstractAlternateCode end
abstract type AbstractTwistedReedSolomonCode <: AbstractLinearCode end

#############################
# concrete types
#############################

#############################
# linearcode.jl
#############################

struct WeightEnumerator
    polynomial::Union{ZZMPolyRingElem,Nemo.AbsSimpleNumFieldElem}
    type::Symbol
end

mutable struct LinearCode <: AbstractLinearCode
    F::CTFieldTypes # base field
    n::Int # length
    k::Int # dimension
    d::Union{Int,Missing} # minimum distance
    l_bound::Int # lower bound on d
    u_bound::Int # upper bound on d
    G::CTMatrixTypes
    H::CTMatrixTypes
    G_stand::CTMatrixTypes
    H_stand::CTMatrixTypes
    P_stand::Union{CTMatrixTypes,Missing} # permutation matrix for G -> G_stand
    weight_enum::Union{WeightEnumerator,Missing}
end

#############################
# MatrixProductCode.jl
#############################

mutable struct MatrixProductCode <: AbstractMatrixProductCode
    F::CTFieldTypes # base field
    n::Int # length
    k::Int # dimension
    d::Union{Int,Missing} # minimum distance
    l_bound::Int # lower bound on d
    u_bound::Int # upper bound on d
    G::CTMatrixTypes
    H::CTMatrixTypes
    G_stand::CTMatrixTypes
    H_stand::CTMatrixTypes
    P_stand::Union{CTMatrixTypes,Missing} # permutation matrix for G -> G_stand
    weight_enum::Union{WeightEnumerator,Missing}
    Cvec::Vector{AbstractLinearCode}
    A::fqPolyRepMatrix
end

#############################
# ReedMuller.jl
#############################

mutable struct ReedMullerCode <: AbstractReedMullerCode
    F::CTFieldTypes
    n::Int # length
    k::Int # dimension
    d::Union{Int,Missing} # minimum distance
    l_bound::Int # lower bound on d
    u_bound::Int # upper bound on d
    r::Integer # order
    m::Integer # number of variables
    G::CTMatrixTypes
    H::CTMatrixTypes
    G_stand::CTMatrixTypes
    H_stand::CTMatrixTypes
    P_stand::Union{CTMatrixTypes,Missing} # permutation matrix for G -> G_stand
    weight_enum::Union{WeightEnumerator,Missing}
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
    d::Union{Int,Missing} # minimum distance
    b::Int # offset
    δ::Int # BCH bound
    HT::Int # Hartmann-Tzeng refinement
    l_bound::Int # lower bound on d
    u_bound::Int # upper bound on d
    qcosets::Vector{Vector{Int}}
    qcosets_reps::Vector{Int}
    def_set::Vector{Int}
    g::CTPolyRingElem
    h::CTPolyRingElem
    e::CTPolyRingElem
    G::CTMatrixTypes
    H::CTMatrixTypes
    G_stand::CTMatrixTypes
    H_stand::CTMatrixTypes
    P_stand::Union{CTMatrixTypes,Missing} # permutation matrix for G -> G_stand
    weight_enum::Union{WeightEnumerator,Missing}
end

mutable struct BCHCode <: AbstractBCHCode
    F::CTFieldTypes # base field
    E::CTFieldTypes # splitting field
    R::CTPolyRing # polynomial ring of generator polynomial
    β::CTFieldElem # n-th root of primitive element of splitting field
    n::Int # length
    k::Int # dimension
    d::Union{Int,Missing} # minimum distance
    b::Int # offset
    δ::Int # BCH bound
    HT::Int # Hartmann-Tzeng refinement
    l_bound::Int # lower bound on d
    u_bound::Int # upper bound on d
    qcosets::Vector{Vector{Int}}
    qcosets_reps::Vector{Int}
    def_set::Vector{Int}
    g::CTPolyRingElem
    h::CTPolyRingElem
    e::CTPolyRingElem
    G::CTMatrixTypes
    H::CTMatrixTypes
    G_stand::CTMatrixTypes
    H_stand::CTMatrixTypes
    P_stand::Union{CTMatrixTypes,Missing} # permutation matrix for G -> G_stand
    weight_enum::Union{WeightEnumerator,Missing}
end

mutable struct ReedSolomonCode <: AbstractReedSolomonCode
    F::CTFieldTypes # base field
    E::CTFieldTypes # splitting field
    R::CTPolyRing # polynomial ring of generator polynomial
    β::CTFieldElem # n-th root of primitive element of splitting field
    n::Int # length
    k::Int # dimension
    d::Union{Int,Missing} # minimum distance
    b::Int # offset
    δ::Int # BCH bound
    HT::Int # Hartmann-Tzeng refinement
    l_bound::Int # lower bound on d
    u_bound::Int # upper bound on d
    qcosets::Vector{Vector{Int}}
    qcosets_reps::Vector{Int}
    def_set::Vector{Int}
    g::CTPolyRingElem
    h::CTPolyRingElem
    e::CTPolyRingElem
    G::CTMatrixTypes
    H::CTMatrixTypes
    G_stand::CTMatrixTypes
    H_stand::CTMatrixTypes
    P_stand::Union{CTMatrixTypes,Missing} # permutation matrix for G -> G_stand
    weight_enum::Union{WeightEnumerator,Missing}
end

#############################
# quasicycliccode.jl
#############################

mutable struct QuasiCyclicCode <: AbstractQuasiCyclicCode
    F::CTFieldTypes # base field
    R::EuclideanRingResidueRing
    n::Int # length
    k::Int # dimension
    d::Union{Int,Missing} # minimum distance
    l_bound::Int # lower bound on d
    u_bound::Int # upper bound on d
    G::Union{CTMatrixTypes,Missing}
    H::Union{CTMatrixTypes,Missing}
    G_stand::Union{CTMatrixTypes,Missing}
    H_stand::Union{CTMatrixTypes,Missing}
    P_stand::Union{CTMatrixTypes,Missing} # permutation matrix for G -> G_stand
    weight_enum::Union{WeightEnumerator,Missing}
    l::Int
    m::Int
    A::MatElem{<:ResElem}
    A_type::Symbol
    W::Matrix{Int}
    type::Int
end

#############################
# GRS_alternate.jl
#############################

mutable struct GeneralizedReedSolomonCode <: AbstractGeneralizedReedSolomonCode
    F::CTFieldTypes # base field
    n::Int # length
    k::Int # dimension
    d::Union{Int,Missing} # minimum distance
    l_bound::Int # lower bound on d
    u_bound::Int # upper bound on d
    scalars::Vector{<:CTFieldElem}
    dual_scalars::Vector{<:CTFieldElem}
    eval_pts::Vector{<:CTFieldElem}
    G::CTMatrixTypes
    H::CTMatrixTypes
    G_stand::CTMatrixTypes
    H_stand::CTMatrixTypes
    P_stand::Union{CTMatrixTypes,Missing} # permutation matrix for G -> G_stand
    weight_enum::Union{WeightEnumerator,Missing} # TODO: should never be missing? is completely known for MDS?
end

mutable struct AlternateCode <: AbstractAlternateCode
    F::CTFieldTypes # base field
    E::CTFieldTypes # extension field
    n::Int # length
    k::Int # dimension
    d::Union{Int,Missing} # minimum distance
    l_bound::Int # lower bound on d
    u_bound::Int # upper bound on d
    scalars::Vector{<:CTFieldElem}
    eval_pts::Vector{<:CTFieldElem}
    G::CTMatrixTypes
    H::CTMatrixTypes
    G_stand::CTMatrixTypes
    H_stand::CTMatrixTypes
    P_stand::Union{CTMatrixTypes,Missing} # permutation matrix for G -> G_stand
    weight_enum::Union{WeightEnumerator,Missing}
end

mutable struct GeneralizedSrivastavaCode <: AbstractGeneralizedSrivastavaCode
    F::CTFieldTypes # base field
    E::CTFieldTypes # extension field
    n::Int # length
    k::Int # dimension
    d::Union{Int,Missing} # minimum distance
    l_bound::Int # lower bound on d
    u_bound::Int # upper bound on d
    a::Vector{<:CTFieldElem}
    w::Vector{<:CTFieldElem}
    z::Vector{<:CTFieldElem}
    t::Int
    G::CTMatrixTypes
    H::CTMatrixTypes
    G_stand::CTMatrixTypes
    H_stand::CTMatrixTypes
    P_stand::Union{CTMatrixTypes,Missing} # permutation matrix for G -> G_stand
    weight_enum::Union{WeightEnumerator,Missing}
end

#############################
# Goppa.jl
#############################

mutable struct GoppaCode <: AbstractGoppaCode
    F::CTFieldTypes # base field
    E::CTFieldTypes # extension field
    n::Int # length
    k::Int # dimension
    d::Union{Int,Missing} # minimum distance
    l_bound::Int # lower bound on d
    u_bound::Int # upper bound on d
    G::CTMatrixTypes
    H::CTMatrixTypes
    G_stand::CTMatrixTypes
    H_stand::CTMatrixTypes
    P_stand::Union{CTMatrixTypes,Missing} # permutation matrix for G -> G_stand
    weight_enum::Union{WeightEnumerator,Missing}
    L::Vector{CTFieldElem}
    g::FqPolyRingElem
end

#############################
# concatenation.jl
#############################

mutable struct ConcatenatedCode <: AbstractLinearCode
    C_in::Union{AbstractLinearCode,Vector{<:AbstractLinearCode}}
    C_out::Union{AbstractLinearCode,Vector{<:AbstractLinearCode}}
    type::Union{Symbol,Vector{Symbol}}
    basis::Union{Missing,Vector{Union{Missing,<:CTFieldElem,Vector{<:CTFieldElem}}}}
    dual_basis::Union{Missing,Vector{Union{Missing,<:CTFieldElem,Vector{<:CTFieldElem}}}}
    F::CTFieldTypes # base field
    n::Int # length
    k::Int # dimension
    d::Union{Int,Missing} # minimum distance
    l_bound::Int # lower bound on d
    u_bound::Int # upper bound on d
    G::CTMatrixTypes
    H::CTMatrixTypes
    G_stand::CTMatrixTypes
    H_stand::CTMatrixTypes
    P_stand::Union{CTMatrixTypes,Missing} # permutation matrix for G -> G_stand
    weight_enum::Union{WeightEnumerator,Missing}
end

#############################
# TwistedReedSolomon.jl
#############################

mutable struct TwistedReedSolomonCode <: AbstractTwistedReedSolomonCode
    F::CTFieldTypes # base field
    n::Int # length
    k::Int # dimension
    d::Union{Int,Missing} # minimum distance
    l_bound::Int # lower bound on d
    u_bound::Int # upper bound on d
    G::CTMatrixTypes
    H::CTMatrixTypes
    G_stand::CTMatrixTypes
    H_stand::CTMatrixTypes
    P_stand::Union{CTMatrixTypes,Missing} # permutation matrix for G -> G_stand
    weight_enum::Union{WeightEnumerator,Missing}
    α::Vector{T} where {T<:CTFieldElem}
    t::Vector{Int}
    h::Vector{Int}
    η::Vector{T} where {T<:CTFieldElem}
    l::Int
end
