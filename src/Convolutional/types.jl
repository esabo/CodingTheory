# Copyright (c) 2024 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
      # abstract types
#############################

abstract type AbstractCode end
abstract type AbstractConvolutionalCode <: AbstractCode end

#############################
      # concrete types
#############################

#############################
   # convolutional_code.jl
#############################

struct PathEnumerator
      polynomial::Union{ZZMPolyRingElem, Nemo.AbsSimpleNumFieldElem}
      type::Symbol
end
    
mutable struct ConvolutionalCode <: AbstractConvolutionalCode
    F::CTFieldTypes # base field
    n::Int # length
    k::Int # dimension
    d::Union{Int, Missing} # free distance
    # l_bound::Int # lower bound on d
    # u_bound::Int # upper bound on d
    D::CTFieldElem # delay operator
    m::Int # memory
    vi::Vector{Int} # constraint lengths
    mnrs::Union{Vector{fqPolyRingElem}, Vector{FqPolyRingElem}}
    int_deg::Int # interal degree
    ext_deg::Int # external degree
    G::CTPolyMatrix
    H::CTPolyMatrix
    path_enum::Union{PathEnumerator, Missing}
end
