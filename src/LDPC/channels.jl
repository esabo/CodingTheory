# Copyright (c) 2023 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

"""
    BinaryErasureChannel(ε::Float64)
    BEC(ε::Float64)

Return the binary erasure channel with erasure probability `ε`.
"""
function BinaryErasureChannel(ε::Float64)
    0 <= ε <= 1 || throw(DomainError("The erasure probability must be in [0, 1]"))
    return BinaryErasureChannel(ε, 1 - ε)
end
BEC(ε::Float64) = BinaryErasureChannel(ε)

"""
    BinarySymmetricChannel(p::Float64)
    BSC(p::Float64)

Return the binary symmetric channel with crossover probability `p`.
"""
function BinarySymmetricChannel(p::Float64)
    0 <= p <= 1 || throw(DomainError("The crossover probability must be in [0, 1]"))
    return BinarySymmetricChannel(p, 1 - _binaryentropy(p))
end
BSC(p::Float64) = BinarySymmetricChannel(p)

"""
    BAWGNChannel(σ::Float64)
    BAWGNC(σ::Float64)

Return the binary (input) additive white Gaussian noise channel with standard deivation `σ`
(noise variance `σ^2`).
"""
function BAWGNChannel(σ::Float64)
    # TODO: is this a good range for this parameter?
    0 <= σ <= 1 || throw(DomainError("The standard deviation must be in [0, 1]"))
    return BAWGNChannel(σ, missing)
end
BAWGNC(σ::Float64) = BAWGNChannel(σ)

# TODO: add an issymmetric parameter to simplify DE later

#############################
      # getter functions
#############################

"""
    erasure_probability(Ch::BinaryErasureChannel)

Return the erasure probability of the binary erasure channel.
"""
erasure_probability(Ch::BinaryErasureChannel) = Ch.param

"""
    crossover_probability(Ch::BinarySymmetricChannel)

Return the crossover probability of the binary symmetric channel.
"""
crossover_probability(Ch::BinarySymmetricChannel) = Ch.param

"""
    standard_deviation(Ch::BAWGNChannel)

Return the standard deviation of the BAWGN channel.
"""
standard_deviation(Ch::BAWGNChannel) = Ch.param

"""
    variance(Ch::BAWGNChannel)

Return the variance of the BAWGN channel.
"""
variance(Ch::BAWGNChannel) = Ch.param^2

#############################
     # general functions
#############################

"""
    capacity(Ch::AbstractClassicalNoiseChannel)

Return the capacity of the noise channel.
"""
function capacity(Ch::AbstractClassicalNoiseChannel)
    ismissing(Ch.capacity) || return Ch.capacity
    
    # TODO: compute capacity functional
    error("Not yet written")
end
