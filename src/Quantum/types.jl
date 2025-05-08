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
abstract type AbstractGeneralizedToricCode <: AbstractStabilizerCodeCSS end

# AbstractQuantumLDPCCode, AbstractQuantumLDPCCSSCode?

abstract type AbstractQuantumNoiseChannel <: AbstractNoiseChannel end

#############################
      # concrete types
#############################

#############################
      # subsystemcode.jl
#############################

# TODO: make sure these have the same info and are in the same order
mutable struct SubsystemCodeCSS <: AbstractSubsystemCodeCSS
      F::CTFieldTypes
      n::Int
      k::Union{Int, Rational{BigInt}}
      r::Int
      d_bare::Union{Int, Missing}
      d_dressed::Union{Int, Missing}
      dx_bare::Union{Int, Missing}
      dx_dressed::Union{Int, Missing}
      dz_bare::Union{Int, Missing}
      dz_dressed::Union{Int, Missing}
      l_bound_bare::Int # lower bound on d_bare
      u_bound_bare::Int # upper bound on d_bare
      l_bound_dressed::Int # lower bound on d_dressed
      u_bound_dressed::Int # upper bound on d_dressed
      l_bound_dx_bare::Int # lower bound on dx_bare
      u_bound_dx_bare::Int # upper bound on dx_bare
      l_bound_dz_bare::Int # lower bound on dz_bare
      u_bound_dz_bare::Int # upper bound on dz_bare
      l_bound_dx_dressed::Int # lower bound on dz_dressed
      u_bound_dx_dressed::Int # upper bound on dz_dressed
      l_bound_dz_dressed::Int # lower bound on dz_dressed
      u_bound_dz_dressed::Int # upper bound on dz_dressed
      stabs::CTMatrixTypes
      X_stabs::CTMatrixTypes
      Z_stabs::CTMatrixTypes
      X_orig_code::Union{S, Missing} where S <: AbstractLinearCode
      Z_orig_code::Union{S, Missing} where S <: AbstractLinearCode
      signs::Vector{zzModRingElem}
      X_signs::Vector{zzModRingElem}
      Z_signs::Vector{zzModRingElem}
      logicals::Vector{Tuple{T, T}} where T <: CTMatrixTypes
      logs_mat::CTMatrixTypes
      char_vec::Vector{zzModRingElem}
      gauge_ops::Vector{Tuple{T, T}} where T <: CTMatrixTypes
      g_ops_mat::CTMatrixTypes
      overcomplete::Bool
      stabs_stand::CTMatrixTypes
      stand_r::Int
      stand_k::Int
      P_stand::Union{CTMatrixTypes, Missing}
      X_metacheck::Union{CTMatrixTypes, Missing}
      Z_metacheck::Union{CTMatrixTypes, Missing}
end

mutable struct SubsystemCode <: AbstractSubsystemCode
      F::CTFieldTypes
      n::Int
      k::Union{Int, Rational{BigInt}}
      r::Int
      d_bare::Union{Int, Missing}
      d_dressed::Union{Int, Missing}
      l_bound_bare::Int # lower bound on d_bare
      u_bound_bare::Int # upper bound on d_bare
      l_bound_dressed::Int # lower bound on d_dressed
      u_bound_dressed::Int # upper bound on d_dressed
      stabs::CTMatrixTypes
      logicals::Vector{Tuple{T, T}} where T <: CTMatrixTypes
      logs_mat::CTMatrixTypes
      char_vec::Vector{zzModRingElem}
      signs::Vector{zzModRingElem}
      gauge_ops::Vector{Tuple{T, T}} where T <: CTMatrixTypes
      g_ops_mat::CTMatrixTypes
      overcomplete::Bool
      stabs_stand::CTMatrixTypes
      stand_r::Int
      stand_k::Int
      P_stand::Union{CTMatrixTypes, Missing}
      metacheck::Union{CTMatrixTypes, Missing}
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
      l_bound::Int # lower bound on d
      u_bound::Int # upper bound on d
      l_bound_dx::Int # lower bound on dx
      u_bound_dx::Int # upper bound on dx
      l_bound_dz::Int # lower bound on dz
      u_bound_dz::Int # upper bound on dz
      stabs::CTMatrixTypes
      X_stabs::CTMatrixTypes
      Z_stabs::CTMatrixTypes
      X_orig_code::Union{S, Missing} where S <: AbstractLinearCode
      Z_orig_code::Union{S, Missing} where S <: AbstractLinearCode
      signs::Vector{zzModRingElem}
      X_signs::Vector{zzModRingElem}
      Z_signs::Vector{zzModRingElem}
      logicals::Vector{Tuple{T, T}} where T <: CTMatrixTypes
      logs_mat::CTMatrixTypes
      char_vec::Vector{zzModRingElem}
      sgn_CWE_stabs::Union{WeightEnumerator, Missing} # signed complete weight enumerator
      sgn_CWE_dual::Union{WeightEnumerator, Missing} # S^⟂
      sgn_CWE_logs::Union{WeightEnumerator, Missing}
      overcomplete::Bool
      pure::Union{Bool, Missing}
      stabs_stand::CTMatrixTypes
      stand_r::Int
      stand_k::Int
      P_stand::Union{CTMatrixTypes, Missing}
      X_metacheck::Union{CTMatrixTypes, Missing}
      Z_metacheck::Union{CTMatrixTypes, Missing}
end

mutable struct StabilizerCode <: AbstractStabilizerCode
      F::CTFieldTypes
      n::Int
      k::Union{Int, Rational{BigInt}}
      d::Union{Int, Missing}
      l_bound::Int # lower bound on d
      u_bound::Int # upper bound on d
      stabs::CTMatrixTypes
      logicals::Vector{Tuple{T, T}} where T <: CTMatrixTypes
      logs_mat::CTMatrixTypes
      char_vec::Vector{zzModRingElem}
      signs::Vector{zzModRingElem}
      sgn_CWE_stabs::Union{WeightEnumerator, Missing} # signed complete weight enumerator
      sgn_CWE_dual::Union{WeightEnumerator, Missing} # S^⟂
      sgn_CWE_logs::Union{WeightEnumerator, Missing}
      overcomplete::Bool
      pure::Union{Bool, Missing}
      stabs_stand::CTMatrixTypes
      stand_r::Int
      stand_k::Int
      P_stand::Union{CTMatrixTypes, Missing}
      metacheck::Union{CTMatrixTypes, Missing}
end

#############################
      # graphstate.jl
#############################

mutable struct GraphStateSubsystem <: AbstractGraphStateSubsystem
      F::CTFieldTypes
      n::Int
      k::Int
      r::Int
      d_bare::Union{Int, Missing}
      d_dressed::Union{Int, Missing}
      l_bound_bare::Int # lower bound on d_bare
      u_bound_bare::Int # upper bound on d_bare
      l_bound_dressed::Int # lower bound on d_dressed
      u_bound_dressed::Int # upper bound on d_dressed
      stabs::CTMatrixTypes
      char_vec::Vector{zzModRingElem}
      signs::Vector{zzModRingElem}
      wtenum::Union{WeightEnumerator, Missing} # signed complete weight enumerator
      overcomplete::Bool
      gauge_ops::Vector{Tuple{T, T}} where T <: CTMatrixTypes
      g_ops_mat::CTMatrixTypes
      stabs_stand::CTMatrixTypes
      stand_r::Int
      stand_k::Int
      P_stand::Union{CTMatrixTypes, Missing}
      metacheck::Union{CTMatrixTypes, Missing}
end

mutable struct GraphStateSubsystemCSS <: AbstractGraphStateSubsystemCSS
      F::CTFieldTypes
      n::Int
      k::Int
      r::Int
      d_bare::Union{Int, Missing}
      d_dressed::Union{Int, Missing}
      dx_bare::Union{Int, Missing}
      dx_dressed::Union{Int, Missing}
      dz_bare::Union{Int, Missing}
      dz_dressed::Union{Int, Missing}
      l_bound_bare::Int # lower bound on d_bare
      u_bound_bare::Int # upper bound on d_bare
      l_bound_dressed::Int # lower bound on d_dressed
      u_bound_dressed::Int # upper bound on d_dressed
      l_bound_dx_bare::Int # lower bound on dx_bare
      u_bound_dx_bare::Int # upper bound on dx_bare
      l_bound_dz_bare::Int # lower bound on dz_bare
      u_bound_dz_bare::Int # upper bound on dz_bare
      l_bound_dx_dressed::Int # lower bound on dz_dressed
      u_bound_dx_dressed::Int # upper bound on dz_dressed
      l_bound_dz_dressed::Int # lower bound on dz_dressed
      u_bound_dz_dressed::Int # upper bound on dz_dressed
      stabs::CTMatrixTypes
      X_stabs::CTMatrixTypes
      Z_stabs::CTMatrixTypes
      X_orig_code::Union{S, Missing} where S <: AbstractLinearCode
      Z_orig_code::Union{S, Missing} where S <: AbstractLinearCode
      signs::Vector{zzModRingElem}
      X_signs::Vector{zzModRingElem}
      Z_signs::Vector{zzModRingElem}
      char_vec::Vector{zzModRingElem}
      sgn_CWE_stabs::Union{WeightEnumerator, Missing} # signed complete weight enumerator
      overcomplete::Bool
      gauge_ops::Vector{Tuple{T, T}} where T <: CTMatrixTypes
      g_ops_mat::CTMatrixTypes
      stabs_stand::CTMatrixTypes
      stand_r::Int
      stand_k::Int
      P_stand::Union{CTMatrixTypes, Missing}
      X_metacheck::Union{CTMatrixTypes, Missing}
      Z_metacheck::Union{CTMatrixTypes, Missing}
end

mutable struct GraphStateStabilizer <: AbstractGraphStateStabilizer
      F::CTFieldTypes
      n::Int
      k::Int
      d::Union{Int, Missing}
      l_bound::Int # lower bound on d
      u_bound::Int # upper bound on d
      stabs::CTMatrixTypes
      char_vec::Vector{zzModRingElem}
      signs::Vector{zzModRingElem}
      sgn_CWE_stabs::Union{WeightEnumerator, Missing} # signed complete weight enumerator
      overcomplete::Bool
      stabs_stand::CTMatrixTypes
      stand_r::Int
      stand_k::Int
      P_stand::Union{CTMatrixTypes, Missing}
      metacheck::Union{CTMatrixTypes, Missing}
end

mutable struct GraphStateStabilizerCSS <: AbstractGraphStateStabilizerCSS
      F::CTFieldTypes
      n::Int
      k::Int
      d::Union{Int, Missing}
      dx::Union{Int, Missing}
      dz::Union{Int, Missing}
      l_bound::Int # lower bound on d
      u_bound::Int # upper bound on d
      l_bound_dx::Int # lower bound on dx
      u_bound_dx::Int # upper bound on dx
      l_bound_dz::Int # lower bound on dz
      u_bound_dz::Int # upper bound on dz
      stabs::CTMatrixTypes
      X_stabs::CTMatrixTypes
      Z_stabs::CTMatrixTypes
      X_orig_code::Union{S, Missing} where S <: AbstractLinearCode
      Z_orig_code::Union{S, Missing} where S <: AbstractLinearCode
      signs::Vector{zzModRingElem}
      X_signs::Vector{zzModRingElem}
      Z_signs::Vector{zzModRingElem}
      char_vec::Vector{zzModRingElem}
      sgn_CWE_stabs::Union{WeightEnumerator,Missing} # signed complete weight enumerator
      overcomplete::Bool
      stabs_stand::CTMatrixTypes
      stand_r::Int
      stand_k::Int
      P_stand::Union{CTMatrixTypes, Missing}
      X_metacheck::Union{CTMatrixTypes, Missing}
      Z_metacheck::Union{CTMatrixTypes, Missing}
end

#############################
# Quantum/product_codes.jl
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
      l_bound::Int # lower bound on d
      u_bound::Int # upper bound on d
      l_bound_dx::Int # lower bound on dx
      u_bound_dx::Int # upper bound on dx
      l_bound_dz::Int # lower bound on dz
      u_bound_dz::Int # upper bound on dz
      stabs::CTMatrixTypes
      X_stabs::CTMatrixTypes
      Z_stabs::CTMatrixTypes
      C1::Union{S, Missing} where S <: AbstractLinearCode
      C2::Union{S, Missing} where S <: AbstractLinearCode
      signs::Vector{zzModRingElem}
      X_signs::Vector{zzModRingElem}
      Z_signs::Vector{zzModRingElem}
      logicals::Vector{Tuple{T, T}} where T <: CTMatrixTypes
      logs_mat::CTMatrixTypes
      char_vec::Vector{zzModRingElem}
      overcomplete::Bool
      stabs_stand::CTMatrixTypes
      stand_r::Int
      stand_k::Int
      P_stand::Union{CTMatrixTypes, Missing}
      sgn_CWE_stabs::Union{WeightEnumerator, Missing} # signed complete weight enumerator
      sgn_CWE_dual::Union{WeightEnumerator, Missing} # S^⟂
      X_metacheck::Union{CTMatrixTypes, Missing}
      Z_metacheck::Union{CTMatrixTypes, Missing}
end

#############################
# Quantum/GeneralizedToricCode.jl
#############################

struct GeneralizedToricCode <: AbstractGeneralizedToricCode
      LR::AbstractAlgebra.Generic.LaurentMPolyWrapRing{fpFieldElem, fpMPolyRing}
      F::CTFieldTypes
      f::CTLRPolyElem
      g::CTLRPolyElem
end

struct FiniteGeneralizedToricCode <: AbstractGeneralizedToricCode
      LR::AbstractAlgebra.Generic.LaurentMPolyWrapRing{fpFieldElem, fpMPolyRing}
      F::CTFieldTypes
      f::CTLRPolyElem
      g::CTLRPolyElem
      a1::Tuple{Int, Int}
      a2::Tuple{Int, Int}
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
