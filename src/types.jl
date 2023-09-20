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
abstract type AbstractConcatenatedCode <: AbstractLinearCode end

abstract type AbstractNoiseChannel end
abstract type AbstractClassicalNoiseChannel <: AbstractNoiseChannel end
abstract type AbstractQuantumNoiseChannel <: AbstractNoiseChannel end
abstract type AbstractBinaryErasureChannel <: AbstractClassicalNoiseChannel end
abstract type AbstractBinarySymmetricChannel <: AbstractClassicalNoiseChannel end
abstract type AbstractBAWGNChannel <: AbstractClassicalNoiseChannel end

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
    l_bound::Int # lower bound on d
    u_bound::Int # upper bound on d
    G::CTMatrixTypes
    H::CTMatrixTypes
    G_stand::CTMatrixTypes
    H_stand::CTMatrixTypes
    P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
    weight_enum::Union{WeightEnumerator, Missing}
end

#############################
       # LDPC/codes.jl
#############################

mutable struct LDPCCode <: AbstractLDPCCode
    F::CTFieldTypes # base field
    n::Int # length
    k::Int # dimension
    d::Union{Int, Missing} # minimum distance
    l_bound::Int # lower bound on d
    u_bound::Int # upper bound on d
    H::CTMatrixTypes
    num_edges::Int
    var_degs::Vector{Int}
    check_degs::Vector{Int}
    col_bound::Int
    row_bound::Int
    limited::Int
    density::Float64
    is_reg::Bool
    Tanner_graph::Union{Figure, Missing}
    λ::fmpq_poly
    ρ::fmpq_poly
    girth::Union{Int, Missing}
    ACEs_per_var_node::Vector{Vector{Int}}
    cycle_lens::Vector{Vector{Int}}
    shortest_cycles::Vector{Vector{Vector{Tuple{Int, Int}}}}
    short_cycle_counts::Dict{Int, Int}
end

#############################
     # LDPC/channels.jl
#############################

struct BinaryErasureChannel <: AbstractBinaryErasureChannel
    param::Float64
    capacity::Float64
end

struct BinarySymmetricChannel <: AbstractBinarySymmetricChannel
    param::Float64
    capacity::Float64
end

mutable struct BAWGNChannel <: AbstractBAWGNChannel
    param::Float64
    capacity::Union{Float64, Missing}
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
    density_evo::Dict{AbstractClassicalNoiseChannel, NTuple{2, Vector{Float64}}}
    threshold::Dict{Type, Float64}
end

#############################
    # MatrixProductCode.jl
#############################

mutable struct MatrixProductCode <: AbstractMatrixProductCode
    F::CTFieldTypes # base field
    n::Int # length
    k::Int # dimension
    d::Union{Int, Missing} # minimum distance
    l_bound::Int # lower bound on d
    u_bound::Int # upper bound on d
    G::CTMatrixTypes
    H::CTMatrixTypes
    G_stand::CTMatrixTypes
    H_stand::CTMatrixTypes
    P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
    weight_enum::Union{WeightEnumerator, Missing}
    Cvec::Vector{AbstractLinearCode}
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
    l_bound::Int # lower bound on d
    u_bound::Int # upper bound on d
    r::Integer # order
    m::Integer # number of variables
    G::CTMatrixTypes
    H::CTMatrixTypes
    G_stand::CTMatrixTypes
    H_stand::CTMatrixTypes
    P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
    weight_enum::Union{WeightEnumerator, Missing}
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
    P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
    weight_enum::Union{WeightEnumerator, Missing}
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
    P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
    weight_enum::Union{WeightEnumerator, Missing}
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
    P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
    weight_enum::Union{WeightEnumerator, Missing}
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
    l_bound::Int # lower bound on d
    u_bound::Int # upper bound on d
    G::Union{CTMatrixTypes, Missing}
    H::Union{CTMatrixTypes, Missing}
    G_stand::Union{CTMatrixTypes, Missing}
    H_stand::Union{CTMatrixTypes, Missing}
    P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
    weight_enum::Union{WeightEnumerator, Missing}
    l::Int
    m::Int
    A::MatElem{<:ResElem}
    A_type::Symbol
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
    l_bound::Int # lower bound on d
    u_bound::Int # upper bound on d
    scalars::Vector{<:CTFieldElem}
    dual_scalars::Vector{<:CTFieldElem}
    eval_pts::Vector{<:CTFieldElem}
    G::CTMatrixTypes
    H::CTMatrixTypes
    G_stand::CTMatrixTypes
    H_stand::CTMatrixTypes
    P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
    weight_enum::Union{WeightEnumerator, Missing} # TODO: should never be missing? is completely known for MDS?
end

#############################
     # concatenation.jl
#############################

mutable struct ConcatenatedCode <: AbstractLinearCode
    C_in::Union{AbstractLinearCode, Vector{<:AbstractLinearCode}}
    C_out::Union{AbstractLinearCode, Vector{<:AbstractLinearCode}}
    type::Union{Symbol, Vector{Symbol}}
    basis::Union{Missing, Vector{Union{Missing, <:CTFieldElem, Vector{<:CTFieldElem}}}}
    dual_basis::Union{Missing, Vector{Union{Missing, <:CTFieldElem, Vector{<:CTFieldElem}}}}
    F::CTFieldTypes # base field
    n::Int # length
    k::Int # dimension
    d::Union{Int, Missing} # minimum distance
    l_bound::Int # lower bound on d
    u_bound::Int # upper bound on d
    G::CTMatrixTypes
    H::CTMatrixTypes
    G_stand::CTMatrixTypes
    H_stand::CTMatrixTypes
    P_stand::Union{CTMatrixTypes, Missing} # permutation matrix for G -> G_stand
    weight_enum::Union{WeightEnumerator, Missing}
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
    X_stabs::CTMatrixTypes
    Z_stabs::CTMatrixTypes
    X_orig_code::Union{LinearCode, Missing}
    Z_orig_code::Union{LinearCode, Missing}
    signs::Vector{nmod}
    X_signs::Vector{nmod}
    Z_signs::Vector{nmod}
    logicals::Vector{Tuple{T, T}} where T <: CTMatrixTypes
    logs_mat::CTMatrixTypes
    char_vec::Vector{nmod}
    gauge_ops::Vector{Tuple{T, T}} where T <: CTMatrixTypes
    g_ops_mat::CTMatrixTypes
    overcomplete::Bool
    stabs_stand::CTMatrixTypes
    stand_r::Int
    stand_k::Int
    P_stand::Union{CTMatrixTypes, Missing}
end
  
mutable struct SubsystemCode <: AbstractSubsystemCode
    F::CTFieldTypes
    n::Int
    k::Union{Int, Rational{BigInt}}
    r::Int
    d::Union{Int, Missing}
    stabs::CTMatrixTypes
    logicals::Vector{Tuple{T, T}} where T <: CTMatrixTypes
    logs_mat::CTMatrixTypes
    char_vec::Vector{nmod}
    signs::Vector{nmod}
    gauge_ops::Vector{Tuple{T, T}} where T <: CTMatrixTypes
    g_ops_mat::CTMatrixTypes
    overcomplete::Bool
    stabs_stand::CTMatrixTypes
    stand_r::Int
    stand_k::Int
    P_stand::Union{CTMatrixTypes, Missing}
end

#############################
      # stabilizercode.jl
#############################

mutable struct StabilizerCodeCSS <: AbstractStabilizerCodeCSS
    F::CTFieldTypes
    n::Int
    k::Union{Int, Rational{BigInt}}
    d::Union{Int, Missing}
    d_x::Union{Int, Missing}
    d_z::Union{Int, Missing}
    stabs::CTMatrixTypes
    X_stabs::CTMatrixTypes
    Z_stabs::CTMatrixTypes
    X_orig_code::Union{LinearCode, Missing}
    Z_orig_code::Union{LinearCode, Missing}
    signs::Vector{nmod}
    X_signs::Vector{nmod}
    Z_signs::Vector{nmod}
    logicals::Vector{Tuple{T, T}} where T <: CTMatrixTypes
    logs_mat::CTMatrixTypes
    char_vec::Vector{nmod}
    sgn_CWE_stabs::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    sgn_CWE_dual::Union{WeightEnumerator, Missing} # S^⟂
    sgn_CWE_logs::Union{WeightEnumerator, Missing}
    overcomplete::Bool
    pure::Union{Bool, Missing}
    stabs_stand::CTMatrixTypes
    stand_r::Int
    stand_k::Int
    P_stand::Union{CTMatrixTypes, Missing}
end
  
mutable struct StabilizerCode <: AbstractStabilizerCode
    F::CTFieldTypes
    n::Int
    k::Union{Int, Rational{BigInt}}
    d::Union{Int, Missing}
    stabs::CTMatrixTypes
    logicals::Vector{Tuple{T, T}} where T <: CTMatrixTypes
    logs_mat::CTMatrixTypes
    char_vec::Vector{nmod}
    signs::Vector{nmod}
    sgn_CWE_stabs::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    sgn_CWE_dual::Union{WeightEnumerator, Missing} # S^⟂
    sgn_CWE_logs::Union{WeightEnumerator, Missing}
    overcomplete::Bool
    pure::Union{Bool, Missing}
    stabs_stand::CTMatrixTypes
    stand_r::Int
    stand_k::Int
    P_stand::Union{CTMatrixTypes, Missing}
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
    char_vec::Vector{nmod}
    signs::Vector{nmod}
    wtenum::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    overcomplete::Bool
    gauge_ops::Vector{Tuple{T, T}} where T <: CTMatrixTypes
    g_ops_mat::CTMatrixTypes
    stabs_stand::CTMatrixTypes
    stand_r::Int
    stand_k::Int
    P_stand::Union{CTMatrixTypes, Missing}
end
  
mutable struct GraphStateSubsystemCSS <: AbstractGraphStateSubsystemCSS
    F::CTFieldTypes
    n::Int
    k::Int
    r::Int
    d::Union{Int, Missing}
    d_x::Union{Int, Missing}
    d_z::Union{Int, Missing}
    stabs::CTMatrixTypes
    X_stabs::CTMatrixTypes
    Z_stabs::CTMatrixTypes
    X_orig_code::Union{LinearCode, Missing}
    Z_orig_code::Union{LinearCode, Missing}
    signs::Vector{nmod}
    X_signs::Vector{nmod}
    Z_signs::Vector{nmod}
    char_vec::Vector{nmod}
    sgn_CWE_stabs::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    overcomplete::Bool
    gauge_ops::Vector{Tuple{T, T}} where T <: CTMatrixTypes
    g_ops_mat::CTMatrixTypes
    stabs_stand::CTMatrixTypes
    stand_r::Int
    stand_k::Int
    P_stand::Union{CTMatrixTypes, Missing}
end
  
mutable struct GraphStateStabilizer <: AbstractGraphStateStabilizer
    F::CTFieldTypes
    n::Int
    k::Int
    d::Union{Int, Missing}
    stabs::CTMatrixTypes
    char_vec::Vector{nmod}
    signs::Vector{nmod}
    sgn_CWE_stabs::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    overcomplete::Bool
    stabs_stand::CTMatrixTypes
    stand_r::Int
    stand_k::Int
    P_stand::Union{CTMatrixTypes, Missing}
end
  
mutable struct GraphStateStabilizerCSS <: AbstractGraphStateStabilizerCSS
    F::CTFieldTypes
    n::Int
    k::Int
    d::Union{Int, Missing}
    d_x::Union{Int, Missing}
    d_z::Union{Int, Missing}
    stabs::CTMatrixTypes
    X_stabs::CTMatrixTypes
    Z_stabs::CTMatrixTypes
    X_orig_code::Union{LinearCode, Missing}
    Z_orig_code::Union{LinearCode, Missing}
    signs::Vector{nmod}
    X_signs::Vector{nmod}
    Z_signs::Vector{nmod}
    char_vec::Vector{nmod}
    sgn_CWE_stabs::Union{WeightEnumerator,Missing} # signed complete weight enumerator
    overcomplete::Bool
    stabs_stand::CTMatrixTypes
    stand_r::Int
    stand_k::Int
    P_stand::Union{CTMatrixTypes, Missing}
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
    d_x::Union{Integer, Missing}
    d_z::Union{Integer, Missing}
    stabs::CTMatrixTypes
    X_stabs::CTMatrixTypes
    Z_stabs::CTMatrixTypes
    C1::Union{LinearCode, Missing}
    C2::Union{LinearCode, Missing}
    signs::Vector{nmod}
    X_signs::Vector{nmod}
    Z_signs::Vector{nmod}
    logicals::Vector{Tuple{T, T}} where T <: CTMatrixTypes
    logs_mat::CTMatrixTypes
    char_vec::Vector{nmod}
    overcomplete::Bool
    stabs_stand::CTMatrixTypes
    stand_r::Int
    stand_k::Int
    P_stand::Union{CTMatrixTypes, Missing}
    sgn_CWE_stabs::Union{WeightEnumerator, Missing} # signed complete weight enumerator
    sgn_CWE_dual::Union{WeightEnumerator, Missing} # S^⟂
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
    hasfield(T, :C) && (C2.C = copy(C.C);)
    if hasfield(T, :Cvec)
        for i in eachindex(C.Cvec)
            C2.C[i] = copy(C.C[i])
        end
    end
    return C2
end
