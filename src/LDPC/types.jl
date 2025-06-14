# Copyright (c) 2023 - 2024 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
# abstract types
#############################

abstract type AbstractLDPCCode <: AbstractLinearCode end

abstract type AbstractNoiseChannel end
abstract type AbstractClassicalNoiseChannel <: AbstractNoiseChannel end
abstract type AbstractBinaryErasureChannel <: AbstractClassicalNoiseChannel end
abstract type AbstractBinarySymmetricChannel <: AbstractClassicalNoiseChannel end
abstract type AbstractBAWGNChannel <: AbstractClassicalNoiseChannel end

#############################
# concrete types
#############################

#############################
# LDPC/codes.jl
#############################

mutable struct LDPCCode <: AbstractLDPCCode
    F::CTFieldTypes # base field
    n::Int # length
    k::Int # dimension
    d::Union{Int,Missing} # minimum distance
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
    # Tanner_graph::Union{Figure, Missing}
    λ::QQPolyRingElem
    ρ::QQPolyRingElem
    girth::Union{Int,Missing}
    ACEs_per_var_node::Vector{Vector{Int}}
    simple_cycles::Vector{Vector{Int}}
    max_cyc_len::Int
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
    capacity::Union{Float64,Missing}
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
    density_evo::Dict{AbstractClassicalNoiseChannel,NTuple{2,Vector{Float64}}}
    threshold::Dict{Type,Float64}
end
