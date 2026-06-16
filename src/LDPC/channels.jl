# Copyright (c) 2023 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

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

function _binary_entropy(p::Float64)
    (p == 0.0 || p == 1.0) && return 0.0
    return -p * log2(p) - (1 - p) * log2(1 - p)
end

"""
$(TYPEDSIGNATURES)

Return the capacity of the noise channel.
"""
capacity(Ch::BinaryErasureChannel) = 1.0 - Ch.ε
capacity(Ch::BinarySymmetricChannel) = 1.0 - _binary_entropy(Ch.p)

"""
$(TYPEDSIGNATURES)

Return the exact BPSK-constrained capacity of the BAWGN channel via numerical integration.
"""
function capacity(Ch::BAWGNChannel)
    σ = Ch.σ
    var = σ^2
    
    # We center the Gaussian integration variable z = y - 1 to make it perfectly 
    # symmetric around 0 for the QuadGK adaptive evaluator.
    function integrand(z)
        # z is the raw Gaussian noise. The received signal is y = 1 + z.
        y = 1.0 + z
        exponent = -2.0 * y / var
        
        # Log2(1 + exp(x)) is numerically unstable for large x.
        # We use a standard numerical trick: 
        # log2(1 + e^x) ≈ x / ln(2) when x > 50
        if exponent > 50.0
            log_term = exponent / log(2.0)
        else
            log_term = log2(1.0 + exp(exponent))
        end
        
        # Multiply by the standard normal PDF
        return exp(-z^2 / (2.0 * var)) / sqrt(2.0 * π * var) * log_term
    end
    
    # Integrate from -∞ to ∞. QuadGK handles infinite bounds via variable substitution natively.
    integral_val, error_bound = quadgk(integrand, -Inf, Inf, rtol=1e-8)
    
    return max(0.0, 1.0 - integral_val)
end

"""
$(TYPEDSIGNATURES)

Return the exact Shannon capacity of the Z-Channel.
"""
function capacity(Ch::ZChannel)
    p = Ch.param
    p == 0.0 && return 1.0
    p == 1.0 && return 0.0
    
    # The analytical maximum mutual information for the Z-channel
    return log2(1.0 + (1.0 - p) * p^(p / (1.0 - p)))
end

"""
$(TYPEDSIGNATURES)

Return the Ergodic Shannon capacity of the Rayleigh Fading channel assuming BPSK modulation.
"""
function capacity(Ch::RayleighFadingChannel)
    # Ergodic capacity: C = ∫ 2a * exp(-a^2) * C_BAWGN(σ / a) da
    # The signal amplitude 'a' effectively scales the noise standard deviation to (σ / a)
    f(a) = 2.0 * a * exp(-a^2) * capacity(BAWGNChannel(Ch.param / a))
    
    # Integrate from 0 to ∞
    cap, _ = quadgk(f, 0.0, Inf; rtol=1e-6)
    return cap
end

"""
$(TYPEDSIGNATURES)

Transmit a binary vector `x` through the channel, returning the corrupted vector `y`.
"""
function transmit(Ch::BinaryErasureChannel, x::AbstractVector{Int})
    # Erasures are conventionally represented by -1 in integer arrays
    return [rand() < Ch.ε ? -1 : bit for bit in x]
end

function transmit(Ch::BinarySymmetricChannel, x::AbstractVector{Int})
    return [rand() < Ch.p ? (1 - bit) : bit for bit in x]
end

function transmit(Ch::BAWGNChannel, x::AbstractVector{Int})
    # BPSK Modulation: 0 -> +1.0, 1 -> -1.0
    bpsk = [bit == 0 ? 1.0 : -1.0 for bit in x]
    noise = randn(length(x)) .* Ch.σ
    return bpsk .+ noise
end

"""
$(TYPEDSIGNATURES)

Simulate transmitting a binary codeword over a Z-Channel.
"""
function transmit(Ch::ZChannel, codeword::Vector{Int})
    p = Ch.param
    # If the bit is 1, it flips to 0 with probability p. 0s stay 0s.
    return [bit == 1 && rand() < p ? 0 : bit for bit in codeword]
end

"""
$(TYPEDSIGNATURES)

Simulate transmitting a binary codeword over a Rayleigh Fading Channel using BPSK.
Returns a tuple `(y, a)` where `y` is the received signal and `a` is the channel state information (fading coefficients).
"""
function transmit(Ch::RayleighFadingChannel, codeword::Vector{Int})
    sigma = Ch.param
    N = length(codeword)
    
    # 1. BPSK Modulation: 0 -> +1.0, 1 -> -1.0
    x = [b == 0 ? 1.0 : -1.0 for b in codeword]
    
    # 2. Generate Rayleigh fading coefficients (normalized so E[a^2] = 1)
    # Using the inverse transform sampling method for perfect analytical Rayleigh generation
    a = sqrt.(.-log.(rand(N)))
    
    # 3. Add AWGN
    noise = sigma .* randn(N)
    
    # Received signal: y = a * x + noise
    y = (a .* x) .+ noise
    
    return y, a
end

"""
$(TYPEDSIGNATURES)

Convert a received channel sequence `y` into Log-Likelihood Ratios (LLRs).
"""
function llr(Ch::BinaryErasureChannel, y::AbstractVector{Int})
    # Erasures (-1) provide 0 information (LLR = 0).
    # Perfect 0 provides +∞, perfect 1 provides -∞.
    inf_val = 1e6 # Cap infinity to prevent numerical overflow in BP
    return [val == -1 ? 0.0 : (val == 0 ? inf_val : -inf_val) for val in y]
end

function llr(Ch::BinarySymmetricChannel, y::AbstractVector{Int})
    # LLR = log((1-p)/p) for a 0, and -log((1-p)/p) for a 1.
    llr_mag = log((1.0 - Ch.p) / Ch.p)
    return [val == 0 ? llr_mag : -llr_mag for val in y]
end

function llr(Ch::BAWGNChannel, y::AbstractVector{Float64})
    # For BPSK over AWGN, the LLR simplifies beautifully to 2y / σ^2
    factor = 2.0 / (Ch.σ^2)
    return y .* factor
end

"""
$(TYPEDSIGNATURES)

Compute the exact Log-Likelihood Ratios (LLRs) for a received Z-Channel signal.
"""
function llr(Ch::ZChannel, received::Vector{Int})
    p = Ch.param
    # If we received a 0:
    # P(y=0|x=0) = 1. P(y=0|x=1) = p. LLR = log(1 / p).
    ll_0 = -log(p) 
    
    # If we received a 1:
    # It is IMPOSSIBLE for x=0 to produce y=1. 
    # Therefore, we are infinitely certain that x=1. LLR = log(0 / (1-p)) = -Inf
    return [r == 0 ? ll_0 : -Inf for r in received]
end

"""
$(TYPEDSIGNATURES)

Compute the exact Log-Likelihood Ratios (LLRs) for a received Rayleigh Fading signal.
Expects `received` to be the tuple `(y, a)` generated by `transmit`.
"""
function llr(Ch::RayleighFadingChannel, received::Tuple{Vector{Float64}, Vector{Float64}})
    y, a = received
    sigma = Ch.param
    
    # For a fading channel with CSI, the LLR formula geometrically scales by 'a'.
    # LLR = 2 * a * y / sigma^2
    multiplier = 2.0 / (sigma^2)
    return multiplier .* a .* y
end

function show(io::IO, Ch::AbstractChannel)
    if isa(Ch, BinaryErasureChannel)
        print(io, "Binary erasure channel with erasure probability $(Ch.param)")
    elseif isa(Ch, BinarySymmetricChannel)
            print(io, "Binary symmetric channel with crossover probability $(Ch.param)")
    elseif isa(Ch, BAWGNChannel)
        print(io, "Binary (input) additive white Gaussian noise channel with standard deviation $(Ch.param)")
    else
        print(io, "Classical noise channel with parameter $(Ch.param)")
    end

    if !ismissing(Ch.capacity)
        println(io, " and capacity $(Ch.capacity).")
    else
        println(io, ".")
    end
end
