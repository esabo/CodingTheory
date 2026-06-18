# Copyright (c) 2023 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # constructors
#############################

# TODO: in tutorial, make sure we explain why we can't pass in L, R here
# TODO: write public conversion functions between all the polynomial types
"""
    LDPCEnsemble(λ::PolyRingElem, ρ::PolyRingElem)

Return the LDPC ensemble determined by the variable degree distribution `λ` and the check
degree distribution `ρ`, both from an edge perspective.
"""
function LDPCEnsemble(λ::PolyRingElem, ρ::PolyRingElem)
    # TODO: check these are proper polynomials, degrees sum to 1, are positive (check if L/R)
    l_avg = _compute_avg_degree(λ)
    r_avg = _compute_avg_degree(ρ)
    L, R = _compute_L_R(λ, ρ, l_avg, r_avg)
    design_rate = Float64(1 - l_avg / r_avg)
    density_evo = Dict{AbstractChannel, NTuple{2, Vector{Float64}}}()
    threshold = Dict{Type, Float64}()
    return LDPCEnsemble(λ, ρ, L, R, Float64(l_avg), Float64(r_avg), design_rate, density_evo,
        threshold)
end

# TODO: ERROR: MethodError: no method matching LDPCEnsemble(::QQPolyRingElem, ::QQPolyRingElem)
"""
    LDPCEnsemble(L::AbstractLDPCCode)

Return the LDPC ensemble determined by the variable degree distribution `λ` and the check
degree distribution `ρ` of `L`, both from an edge perspective.
"""
LDPCEnsemble(L::AbstractLDPCCode) = LDPCEnsemble(L.λ, L.ρ)

#############################
      # getter functions
#############################

# do we not have getter functions for this type?

#############################
      # setter functions
#############################

#############################
     # general functions
#############################

# TODO: are _d_poly(...) and _integrate_poly(...) useful? currently unused, probably delete
_d_poly(vec::Vector{<:Real}) = [vec[i] * (i - 1) for i in 2:length(vec)]
_d_poly(f::PolyRingElem) = derivative(f)

_integrate_poly(vec::Vector{T}) where T <: Real = [zero(T); [c / i for (i, c) in enumerate(vec)]]
_integrate_poly(f::PolyRingElem) = integral(f)

# _poly_eval, _d_poly_eval, and _integrate_poly_0_1 are all used
_poly_eval(x::Real, vec::Vector{<:Real}) = sum(c * x^(i - 1) for (i, c) in enumerate(vec))
_poly_eval(x::Real, f::PolyRingElem) = _poly_eval(x, Float64.(coeff.(f, 0:degree(f))))

_d_poly_eval(x::Real, vec::Vector{<:Real}) = sum((i - 1) * c * x^(i - 2) for (i, c) in enumerate(vec))
_d_poly_eval(x::Real, f::PolyRingElem) = _d_poly_eval(x, Float64.(coeff.(f, 0:degree(f))))

_integrate_poly_0_1(vec::Vector{<:Real}) = sum(c / i for (i, c) in enumerate(vec))
_integrate_poly_0_1(f::PolyRingElem) = _integrate_poly_0_1(Float64.(coeff.(f, 0:degree(f))))

# TODO: move to utils, check if already there, export
# possibly just go ahead and extend to the nonbinary case then call with a 2 here
_binary_entropy(x::Real) = x * log2(1 / x) + (1 - x) * log2(1 / (1 - x))

# these should all be useful. Note that the QQ version of _compute_λ_ρ is necessary for the output to also be a QQPoly. The non-QQ version is for when the poly isn't directly callable (as in RealPoly)
_compute_avg_degree(f::PolyRingElem) = inv(sum(coeff(f, i) / (i + 1) for i in 0:degree(f)))
_compute_L_R(λ::PolyRingElem, ρ::PolyRingElem, l_avg, r_avg) = (integral(λ) * l_avg, integral(ρ) * r_avg)
_compute_L_R(λ::PolyRingElem, ρ::PolyRingElem) = _compute_L_R(λ, ρ, _compute_avg_degree(λ), _compute_avg_degree(ρ))
_compute_λ_ρ(L::PolyRingElem, R::PolyRingElem) = (derivative(L) / _d_poly_eval(1, L), derivative(R) / _d_poly_eval(1, R))
_compute_λ_ρ(L::QQPolyRingElem, R::QQPolyRingElem) = (derivative(L) / derivative(L)(1), derivative(R) / derivative(R)(1))

function _L2_dist_sq(p1::Vector{Float64}, p2::Vector{Float64})
    @assert length(p1) == length(p2)
    v = p1 .- p2
    v2 = [sum(v[j] * v[k + 1 - j] for j in max(1, k + 1 - length(v)):min(k, length(v))) for k in 1:2length(v) - 1]
    return _integrate_poly_0_1(v2)
end

# Internal helper to grab the defining parameter (ε, p, σ, etc.)
_param(Ch::AbstractChannel) = getfield(Ch, 1)

# Include the h::UInt salt for proper Julia hashing performance
Base.hash(Ch::AbstractChannel, h::UInt) = hash(_param(Ch), hash(typeof(Ch), h))

# Type-stable equality check: channels are only equal if they are the exact same type AND have the same parameter
Base.isequal(Ch1::T, Ch2::T) where {T <: AbstractChannel} = isequal(_param(Ch1), _param(Ch2))
Base.isequal(::AbstractChannel, ::AbstractChannel) = false

# It is also best practice in Julia to map `==` to `isequal` for custom types
Base.:(==)(Ch1::AbstractChannel, Ch2::AbstractChannel) = isequal(Ch1, Ch2)

# function Base.setproperty!(Ch::BAWGNChannel, key, val)
#     key == :capacity && (setfield!(Ch, key, val);)
#     key == :capacity || @warn "Channel not updated. Create a new channel instead of changing the noise on an existing channel."
# end

function _density_evolution!(E::LDPCEnsemble, Ch::AbstractChannel)
    if isa(Ch, BinaryErasureChannel)
        λ_vec = Float64.(coeff.(E.λ, 0:degree(E.λ)))
        ρ_vec = Float64.(coeff.(E.ρ, 0:degree(E.ρ)))
        E.density_evo[Ch] = _density_evolution_BEC(λ_vec, ρ_vec, Ch.ε)
    else
        error("Only BEC has been implemented so far")
    end
    return nothing
end

"""
    density_evolution(E::LDPCEnsemble, Ch::AbstractChannel)

Return the density evolution of the LDPC ensemble given the noise channel.
"""
function density_evolution(E::LDPCEnsemble, Ch::AbstractChannel)
    Ch ∈ keys(E.density_evo) || _density_evolution!(E::LDPCEnsemble, Ch::AbstractChannel)
    return E.density_evo[Ch]
end

function _density_evolution_BEC(λ::Vector{<:Real}, ρ::Vector{<:Real}, ε::Real;
    max_iters::Int=500, tol::Float64=1e-9)

    iter = 0
    evo_x = [ε]
    evo_y = [1.0]
    while evo_x[end] > tol && iter < max_iters
        iter += 1
        push!(evo_y, 1 - _poly_eval(1 - evo_x[end], ρ))
        push!(evo_x, ε * _poly_eval(evo_y[end], λ))
    end
    return evo_x, evo_y
end

"""
$(TYPEDSIGNATURES)

Return the multiplicative gap of the ensemble with respect to the given channel.
"""
function multiplicative_gap(E::LDPCEnsemble, Ch::AbstractChannel)
    if !haskey(E.threshold, typeof(Ch))
        # For BEC, we can use the exact analytical threshold.
        if typeof(Ch) == BinaryErasureChannel
            # FIXED: Pass the ensemble E and the type of the channel
            E.threshold[typeof(Ch)] = optimal_threshold(E, typeof(Ch))
        else
            error("Threshold computation for this channel is not yet implemented.")
        end
    end
    
    thresh = E.threshold[typeof(Ch)]
    C_thresh = capacity(typeof(Ch)(thresh))
    
    return (C_thresh - E.design_rate) / C_thresh
end

"""
    multiplicative_gap_lower_bound(E::LDPCEnsemble)

Return a lower bound on the multiplicative gap of the ensemble
"""
multiplicative_gap_lower_bound(E::LDPCEnsemble) = (E.design_rate^E.r_avg * (1 - E.design_rate)) / (1 + E.design_rate^E.r_avg * (1 - E.design_rate))

"""
    density_lower_bound(Ch::AbstractChannel, gap::Real)

Return a lower bound on the density of a (full rank) parity-check matrix for the channel
given the multiplicative gap.
"""
function density_lower_bound(Ch::AbstractChannel, gap::Real)
    0 < gap < 1 || throw(DomainError("Multiplicative gap should be in (0, 1)"))
    if isa(Ch, BinaryErasureChannel)
        temp = log(1 - Ch.ε)
        return (Ch.ε * (log(gap) - (log(Ch.ε) - temp))) / ((1 - Ch.ε) * (1 - gap) * temp)
    else
        @error "Not yet implemented"
    end
end

"""
    check_concentrated_degree_distribution(Ch::BinaryErasureChannel, gap::Real)

Return the check-concentrated degree distribution `(λ, ρ)` for the binary erasure channel
given the desired multiplicative gap.
"""
function check_concentrated_degree_distribution(Ch::BinaryErasureChannel, gap::Real)
    # Euler-Mascheroni constant
    γ = 0.577215664901533
    pi26 = π^2 / 6
    temp = 1 - Ch.param
    c = temp^pi26 * exp((pi26 - γ) * Ch.param)
    N = max(ceil(Int, 1 - c * temp * (1 - gap) / gap), ceil(Int, temp^(-1 / Ch.param)))
    α = log(1 / temp) / log(N)

    # this check isn't going to work in most cases
    # isinteger(1 / α) || error("1/α is not an integer")

    # this fixed choice does give the correct answer (p 115)
    # N = 13
    # α = 1/5

    λ_vec = zeros(N - 1)
    λ_vec[1] = α
    for i in 2:N - 1
        λ_vec[i] = ((i - 1) / i) * (1 - α / (i - 1)) * λ_vec[i - 1]
    end
    norm = sum(λ_vec)
    λ_vec ./= norm

    _, x = PolynomialRing(RealField(), :x)
    λ = sum(λi * x^i for (i, λi) in enumerate(λ_vec))
    ρ = x^round(Int, 1 / α)
    return λ, ρ
end

"""
    optimal_lambda(ρ, l_max, param, var_type; Δ = 1e-3)

Find the optimal variable node distribution `λ` given the check node
distribution `ρ`, maximum variable node degree `l_max`, and target parameter
`param` which refers to threshold if `var_type == :ε` or rate if `var_type == :r`.

# Notes
* `Δ` refers to the step size for `x` in the LP to solve for `λ`.
"""
function optimal_lambda end

"""
    optimal_rho(λ, r_max, param, var_type; Δ = 1e-3)

Find the optimal check node distribution `ρ` given the variable node
distribution `λ`, maximum check node degree `r_max`, and target parameter `param`
which refers to threshold if `var_type == :ε` or rate if `var_type == :r`.

# Notes
* `Δ` refers to the step size for x in the LP to solve for ρ.
"""
function optimal_rho end

"""
    optimal_lambda_and_rho(l_max, r_max, param, var_type; Δρ = 1e-2, Δλ = 1e-3)

Find the optimal distribution pair λ, ρ given the `param`, where `param` is
either a threshold if `var_type == :ε` or a target rate if `var_type == :r`.

# Notes
* `Δρ` gives the step size for possible values of `c` where
    ρ = (1 - c) * x^(r_max - 2) + c * x^(r_max - 1)
* `Δλ` gives the step size for values of `x` in the LP for finding λ given ρ.
"""
function optimal_lambda_and_rho end

# Fast analytical approximation of the GA check-node evolution function (Chung 2001)
function _phi_GA(x::Float64)
    x <= 0.0 && return 1.0
    x > 10.0 && return sqrt(pi / x) * exp(-x / 4.0) * (1.0 - 10.0 / (7.0 * x))
    return exp(-0.4527 * x^0.86 + 0.0218)
end

# Inverse of the phi function via fast binary search
function _inv_phi_GA(y::Float64)
    y >= 1.0 && return 0.0
    y <= 0.0 && return typemax(Float64)
    
    # Binary search bounds
    low = 0.0
    high = 100.0 # High enough to act as infinity for typical LLR means
    
    # Expand high bound if necessary
    while _phi_GA(high) > y
        high *= 2.0
    end
    
    # Fast bisection
    for _ in 1:50
        mid = (low + high) / 2.0
        if _phi_GA(mid) > y
            low = mid
        else
            high = mid
        end
    end
    return (low + high) / 2.0
end

"""
    _density_evolution_GA(λ_vec, ρ_vec, σ; max_iters=500, tol=1e-6)

Perform Gaussian Approximation Density Evolution for the BAWGN channel.
Returns `true` if the DE successfully decodes (mean LLR approaches infinity), `false` otherwise.
"""
function _density_evolution_GA(λ_vec::Vector{<:Real}, ρ_vec::Vector{<:Real}, σ::Float64; max_iters::Int=500, tol::Float64=1e-6)
    # Initial channel LLR mean for BPSK over AWGN: 2 / σ^2
    m_u0 = 2.0 / (σ^2)
    m_v = m_u0
    
    for iter in 1:max_iters
        # 1. Check Node Update: Expected mean from check to variable
        # sum_{j} ρ_j * phi_inv( 1 - [1 - phi(m_v)]^(j-1) )
        term_phi = _phi_GA(m_v)
        m_u = 0.0
        for (j, rho_j) in enumerate(ρ_vec)
            if rho_j > 0
                power_val = (1.0 - term_phi)^(j - 1)
                m_u += rho_j * _inv_phi_GA(1.0 - power_val)
            end
        end
        
        # 2. Variable Node Update: Expected mean from variable to check
        # m_v = m_u0 + sum_{i} λ_i * (i-1) * m_u
        m_v_next = m_u0
        for (i, lam_i) in enumerate(λ_vec)
            if lam_i > 0
                m_v_next += lam_i * (i - 1) * m_u
            end
        end
        
        # Check for convergence (if mean LLR is growing massively, it decoded)
        if m_v_next > 50.0 
            return true
        end
        
        # Check if it stalled (error floor / waterfall failure)
        if abs(m_v_next - m_v) < tol
            return false
        end
        
        m_v = m_v_next
    end
    return false
end

"""
$(TYPEDSIGNATURES)

Return the optimal threshold for the LDPC ensemble over the specified Channel type.
"""
function optimal_threshold(E::LDPCEnsemble, ::Type{BinaryErasureChannel}; Δ::Float64=1e-4)
    xs = Δ:Δ:1.0
    return minimum(x / (_poly_eval(1.0 - _poly_eval(1.0 - x, E.ρ), E.λ)) for x in xs)
end

function optimal_threshold(E::LDPCEnsemble, ::Type{BAWGNChannel}; tol::Float64=1e-4)
    λ_vec = Float64.(coeff.(E.λ, 0:degree(E.λ)))
    ρ_vec = Float64.(coeff.(E.ρ, 0:degree(E.ρ)))
    
    # Binary search for the maximum noise standard deviation (σ) that still decodes
    # Note: Higher σ is a WORSE channel (unlike BEC where higher ε is worse)
    low_σ = 0.1  # Very good channel (should decode)
    high_σ = 3.0 # Very bad channel (should fail)
    
    # Ensure our bounds are valid
    while _density_evolution_GA(λ_vec, ρ_vec, high_σ)
        high_σ *= 2.0
    end
    
    # Fast bisection to find the threshold
    for _ in 1:50
        mid_σ = (low_σ + high_σ) / 2.0
        if _density_evolution_GA(λ_vec, ρ_vec, mid_σ)
            low_σ = mid_σ # Decoded! Can we handle more noise?
        else
            high_σ = mid_σ # Failed! Need less noise.
        end
        
        if (high_σ - low_σ) < tol
            break
        end
    end
    
    return low_σ
end

"""
$(TYPEDSIGNATURES)

Generate the `(x, y)` curve data for the EXIT chart of the ensemble over a given channel.
Returns a tuple `(vnd_x, vnd_y, cnd_x, cnd_y)`.

# Notes
* `vnd` curves represent the Variable Node Decoder mutual information transfer.
* `cnd` curves represent the Check Node Decoder mutual information transfer.
* On a standard EXIT chart, the axes are swapped for the CND curve to visualize the 
  decoding tunnel.
"""
function EXIT_chart_data(E::LDPCEnsemble, Ch::BinaryErasureChannel; pts::Int=100)
    # Mutual Information always sweeps from 0.0 to 1.0
    I_A = collect(range(0.0, 1.0, length=pts))
    
    λ_vec = Float64.(coeff.(E.λ, 0:degree(E.λ)))
    ρ_vec = Float64.(coeff.(E.ρ, 0:degree(E.ρ)))
    ε = Ch.ε
    
    # 1. Variable Node Curve: I_E = 1 - ε * λ(1 - I_A)
    # Plotted normally: x = I_A, y = I_E
    vnd_x = I_A
    vnd_y = [1.0 - ε * _poly_eval(1.0 - ia, λ_vec) for ia in I_A]
    
    # 2. Check Node Curve: I_E = ρ(I_A)
    # Plotted inverted to form the tunnel: x = I_E, y = I_A
    cnd_x = [_poly_eval(ia, ρ_vec) for ia in I_A]
    cnd_y = I_A
    
    return vnd_x, vnd_y, cnd_x, cnd_y
end

"""
    EXIT_chart_plot(E::LDPCEnsemble, Ch::AbstractChannel; tol::Float64 = 1e-9)

Return a plot of the EXIT chart for the ensemble given the channel up to a numerical tolerance of `tol`.

# Note
- Run `using Makie` to activate this extension.
"""
function EXIT_chart_plot end

# Map our existing Gaussian Approximation functions to Mutual Information
_J_MI(σ::Float64) = 1.0 - _phi_GA((σ^2) / 2.0)
_inv_J_MI(I::Float64) = sqrt(2.0 * _inv_phi_GA(clamp(1.0 - I, 0.0, 1.0)))

"""
$(TYPEDSIGNATURES)

Run a Protograph EXIT (PEXIT) analysis on the base matrix `B` for a given channel 
LLR standard deviation `sigma_ch`. 

# Arguments
* `B::Matrix{Int}`: The protograph base matrix.
* `sigma_ch::Vector{Float64}`: The channel LLR standard deviation for each variable node. 
  If a node is punctured, its value should be strictly `0.0`.

# Returns
* `(decoded::Bool, I_APP::Vector{Float64}, iters::Int)`
"""
function _PEXIT_AWGN(B::Matrix{Int}, sigma_ch::Vector{Float64}; max_iters::Int=1000, tol::Float64=1e-4)
    m, n = size(B)
    length(sigma_ch) == n || throw(ArgumentError("sigma_ch must match the number of variable nodes in B"))
    
    # State matrices for the Extrinsic Mutual Information tracking
    I_V2C = zeros(Float64, m, n)
    I_C2V = zeros(Float64, m, n)
    I_APP = zeros(Float64, n)
    
    # Initialize Variable-to-Check messages purely with channel info
    for j in 1:n
        if sigma_ch[j] > 0.0
            mi_ch = _J_MI(sigma_ch[j])
            for i in 1:m
                if B[i,j] > 0
                    I_V2C[i,j] = mi_ch
                end
            end
        end
    end
    
    for iter in 1:max_iters
        # 1. Check Node Update (combines MI in the dual domain)
        for i in 1:m
            for j in 1:n
                B[i,j] == 0 && continue
                
                sum_inv_J = 0.0
                for k in 1:n
                    if B[i,k] > 0
                        # Maintain extrinsic principle: exclude one edge connecting to j
                        edges = (k == j) ? B[i,k] - 1 : B[i,k]
                        if edges > 0
                            sum_inv_J += edges * (_inv_J_MI(1.0 - I_V2C[i,k]))^2
                        end
                    end
                end
                I_C2V[i,j] = 1.0 - _J_MI(sqrt(sum_inv_J))
            end
        end
        
        # 2. Variable Node Update (combines MI in the variance domain)
        for j in 1:n
            for i in 1:m
                B[i,j] == 0 && continue
                
                sum_inv_J = sigma_ch[j]^2
                for k in 1:m
                    if B[k,j] > 0
                        # Maintain extrinsic principle: exclude one edge connecting to i
                        edges = (k == i) ? B[k,j] - 1 : B[k,j]
                        if edges > 0
                            sum_inv_J += edges * (_inv_J_MI(I_C2V[k,j]))^2
                        end
                    end
                end
                I_V2C[i,j] = _J_MI(sqrt(sum_inv_J))
            end
        end
        
        # 3. Calculate A Posteriori Probability (APP) MI for convergence check
        all_decoded = true
        max_delta = 0.0
        
        for j in 1:n
            sum_inv_J = sigma_ch[j]^2
            for k in 1:m
                if B[k,j] > 0
                    sum_inv_J += B[k,j] * (_inv_J_MI(I_C2V[k,j]))^2
                end
            end
            
            new_I_APP = _J_MI(sqrt(sum_inv_J))
            max_delta = max(max_delta, abs(new_I_APP - I_APP[j]))
            I_APP[j] = new_I_APP
            
            # Punctured nodes rely entirely on graph edges, so we ensure ALL nodes reach 1.0
            if I_APP[j] < 0.999
                all_decoded = false
            end
        end
        
        # Early stopping if completely decoded, or if stalled (error floor hit)
        if all_decoded
            return true, I_APP, iter
        end
        if max_delta < tol
            return false, I_APP, iter
        end
    end
    
    return false, I_APP, max_iters
end

"""
$(TYPEDSIGNATURES)

Find the exact AWGN noise threshold (maximum standard deviation `σ_n`) for a 
protograph base matrix using PEXIT analysis.

# Arguments
* `B::Matrix{Int}`: The protograph base matrix.
* `punctured::Vector{Bool}`: A boolean vector indicating if a column is punctured. 
  Defaults to all `false`.
"""
function protograph_threshold(B::Matrix{Int}; punctured::Union{Vector{Bool}, Nothing} = nothing)
    m, n = size(B)
    if isnothing(punctured)
        punctured = fill(false, n)
    end
    length(punctured) == n || throw(ArgumentError("punctured array length must match number of columns in B"))
    
    # Binary search bounds for the channel noise standard deviation (σ_n)
    low_σ_n = 0.1  # Very low noise -> Should trivially decode
    high_σ_n = 3.0 # Very high noise -> Should fail
    
    # Helper to generate the LLR standard deviation vector. 
    # For AWGN, LLR variance = 4 / σ_n^2, so LLR standard deviation = 2 / σ_n
    function _get_sigma_ch(σ_n::Float64)
        sig_ch = zeros(Float64, n)
        for j in 1:n
            if !punctured[j]
                sig_ch[j] = 2.0 / σ_n
            end
        end
        return sig_ch
    end
    
    # Ensure our upper bound is actually failing
    while _PEXIT_AWGN(B, _get_sigma_ch(high_σ_n))[1]
        high_σ_n *= 2.0
    end
    
    # Fast bisection
    for _ in 1:50
        mid_σ_n = (low_σ_n + high_σ_n) / 2.0
        sig_ch = _get_sigma_ch(mid_σ_n)
        
        decoded, _, _ = _PEXIT_AWGN(B, sig_ch)
        if decoded
            low_σ_n = mid_σ_n  # Matrix survived this noise, can we push it harder?
        else
            high_σ_n = mid_σ_n # Matrix failed, lower the noise
        end
        
        if (high_σ_n - low_σ_n) < 1e-4
            break
        end
    end
    
    return low_σ_n
end

"""
$(TYPEDSIGNATURES)

Generate the `(x, y)` curve data for the averaged EXIT chart of a Protograph base matrix `B`.
Returns a tuple `(vnd_x, vnd_y, cnd_x, cnd_y)` which seamlessly plugs into standard plotting functions.

# Arguments
* `B::Matrix{Int}`: The protograph base matrix.
* `sigma_ch::Vector{Float64}`: The channel LLR standard deviation for each variable node.
"""
function PEXIT_chart_data(B::Matrix{Int}, sigma_ch::Vector{Float64}; pts::Int=100)
    m, n = size(B)
    total_edges = sum(B)
    
    I_A_sweep = collect(range(0.0, 1.0, length=pts))
    
    vnd_x = I_A_sweep
    vnd_y = zeros(Float64, pts)
    cnd_x = zeros(Float64, pts)
    cnd_y = I_A_sweep
    
    for (idx, I_A) in enumerate(I_A_sweep)
        # 1. Variable Node Curve (VND) Average Transfer
        sum_I_E_vnd = 0.0
        for j in 1:n
            # Total edges connected to variable node j
            deg_v = sum(B[:, j]) 
            for i in 1:m
                B[i,j] == 0 && continue
                
                # Variance coming in: Channel + (All other edges from checks) * inv_J(I_A)^2
                edges_from_other_checks = deg_v - 1 
                var_in = sigma_ch[j]^2 + edges_from_other_checks * (_inv_J_MI(I_A))^2
                
                I_E = _J_MI(sqrt(var_in))
                sum_I_E_vnd += B[i,j] * I_E # Weight by the number of edges
            end
        end
        vnd_y[idx] = sum_I_E_vnd / total_edges
        
        # 2. Check Node Curve (CND) Average Transfer
        sum_I_E_cnd = 0.0
        for i in 1:m
            # Total edges connected to check node i
            deg_c = sum(B[i, :])
            for j in 1:n
                B[i,j] == 0 && continue
                
                edges_from_other_vars = deg_c - 1
                var_in = edges_from_other_vars * (_inv_J_MI(1.0 - I_A))^2
                
                I_E = 1.0 - _J_MI(sqrt(var_in))
                sum_I_E_cnd += B[i,j] * I_E
            end
        end
        cnd_x[idx] = sum_I_E_cnd / total_edges
    end
    
    return vnd_x, vnd_y, cnd_x, cnd_y
end

"""
$(TYPEDSIGNATURES)

Generate the `(x, y)` curve data for the EXIT chart of the ensemble over a BAWGN channel
using Gaussian Approximation.
"""
function EXIT_chart_data(E::LDPCEnsemble, Ch::BAWGNChannel; pts::Int=100)
    # For AWGN, LLR variance = 4 / σ_n^2. Thus, LLR std dev = 2 / σ_n
    sigma_ch = 2.0 / Ch.σ # FIXED: Use Ch.σ instead of Ch.param
    return _EXIT_chart_GA(E, sigma_ch, pts)
end

"""
$(TYPEDSIGNATURES)

Generate the `(x, y)` curve data for the EXIT chart of the ensemble over any arbitrary 
symmetric channel using the AWGN-equivalent capacity approximation.
"""
function EXIT_chart_data(E::LDPCEnsemble, Ch::AbstractChannel; pts::Int=100)
    # 1. Find the exact Mutual Information (Capacity) of the channel
    I_ch = capacity(Ch)
    
    # 2. Map it to an AWGN-equivalent LLR standard deviation
    sigma_ch = _inv_J_MI(I_ch)
    
    return _EXIT_chart_GA(E, sigma_ch, pts)
end

"""
$(TYPEDSIGNATURES)

Internal engine to compute the GA EXIT chart curves given an initial channel LLR standard deviation.
"""
function _EXIT_chart_GA(E::LDPCEnsemble, sigma_ch::Float64, pts::Int)
    I_A_sweep = collect(range(0.0, 1.0, length=pts))
    
    λ_vec = Float64.(coeff.(E.λ, 0:degree(E.λ)))
    ρ_vec = Float64.(coeff.(E.ρ, 0:degree(E.ρ)))
    
    vnd_x = I_A_sweep
    vnd_y = zeros(Float64, pts)
    cnd_x = zeros(Float64, pts)
    cnd_y = I_A_sweep
    
    for (idx, I_A) in enumerate(I_A_sweep)
        # --- Variable Node Curve (VND) ---
        # I_E = sum( λ_i * J( sqrt( sigma_ch^2 + (i-1) * inv_J(I_A)^2 ) ) )
        inv_J_A_sq = (_inv_J_MI(I_A))^2
        sum_vnd = 0.0
        
        for (i, lam_i) in enumerate(λ_vec)
            if lam_i > 0
                var_in = sigma_ch^2 + (i - 1) * inv_J_A_sq
                sum_vnd += lam_i * _J_MI(sqrt(var_in))
            end
        end
        vnd_y[idx] = sum_vnd
        
        # --- Check Node Curve (CND) ---
        # I_E = 1 - sum( ρ_j * J( sqrt( (j-1) * inv_J(1 - I_A)^2 ) ) )
        inv_J_1_minus_A_sq = (_inv_J_MI(1.0 - I_A))^2
        sum_cnd = 0.0
        
        for (j, rho_j) in enumerate(ρ_vec)
            if rho_j > 0
                var_in = (j - 1) * inv_J_1_minus_A_sq
                sum_cnd += rho_j * _J_MI(sqrt(var_in))
            end
        end
        cnd_x[idx] = 1.0 - sum_cnd
    end
    
    return vnd_x, vnd_y, cnd_x, cnd_y
end

# Standard Q-function: Tail probability of the standard normal distribution
_Q_function(x::Real) = 0.5 * erfc(x / sqrt(2.0))

"""
$(TYPEDSIGNATURES)

Estimate the Block Error Rate (BLER) of an LDPC code at a finite block length `n` 
using the refined scaling law (Amraoui et al.).

# Arguments
* `threshold`: The theoretical asymptotic decoding limit (e.g., `\\epsilon^*` for BEC or `\\sigma^*` for AWGN).
* `param`: The actual channel parameter currently being evaluated.
* `n`: The physical block length of the codeword.
* `alpha`: The scaling parameter controlling the waterfall slope (variance of the decoding trajectory).
* `beta`: The shift parameter controlling the finite-length performance penalty.

# Notes
* The effective finite-length threshold is mathematically shifted by `beta * n^(-2/3)`.
"""
function finite_length_estimate(threshold::Float64, param::Float64, n::Int, alpha::Float64, beta::Float64)
    # 1. Calculate the effective threshold at this specific block length
    # The term n^(-2/3) dictates the exact shift penalty.
    effective_threshold = threshold - beta * (n ^ (-2.0/3.0))
    
    # 2. Calculate the distance from the new effective threshold
    # Positive delta means we are in the "good" channel region.
    delta = effective_threshold - param
    
    # 3. Scale by the slope parameter alpha and sqrt(n)
    z = (sqrt(n) * delta) / alpha
    
    return _Q_function(z)
end

"""
$(TYPEDSIGNATURES)

Run Multi-Edge Type (MET) Density Evolution using Gaussian Approximation for an AWGN channel.
Returns `true` if the ensemble decodes, `false` otherwise.

# Arguments
* `E::METEnsemble`: The Multi-Edge Type ensemble definition.
* `sigma_ch::Float64`: The AWGN channel LLR standard deviation (2.0 / sigma_noise).
"""
function density_evolution_MET_GA(E::METEnsemble, sigma_ch::Float64; max_iters::Int=1000, tol::Float64=1e-5)
    ne = E.num_edge_types
    nv_classes = size(E.var_profiles, 1)
    nc_classes = size(E.chk_profiles, 1)
    
    # I_V2C[e] = The average Mutual Information traveling along edge type 'e' from Var to Check
    I_V2C = zeros(Float64, ne)
    # I_C2V[e] = The average Mutual Information traveling along edge type 'e' from Check to Var
    I_C2V = zeros(Float64, ne)
    
    # Initialize V2C messages with purely channel information
    for e in 1:ne
        var_sum = 0.0
        weight_sum = 0.0
        for v in 1:nv_classes
            frac = E.var_profiles[v, 1]
            is_transmitted = E.var_profiles[v, 2]
            edges_of_type_e = E.var_profiles[v, 2 + e]
            
            if edges_of_type_e > 0
                # If transmitted, it has channel variance. If punctured (0.0), it has 0 channel variance.
                ch_var = (is_transmitted * sigma_ch)^2
                var_sum += frac * edges_of_type_e * ch_var
                weight_sum += frac * edges_of_type_e
            end
        end
        # Map the average starting variance to Mutual Information
        if weight_sum > 0
            I_V2C[e] = _J_MI(sqrt(var_sum / weight_sum))
        end
    end
    
    for iter in 1:max_iters
        # ---------------------------------------------------------
        # 1. Check Node Update (CND): Combine MI in the dual domain
        # ---------------------------------------------------------
        new_I_C2V = zeros(Float64, ne)
        for e in 1:ne
            mi_sum = 0.0
            weight_sum = 0.0
            
            for c in 1:nc_classes
                frac = E.chk_profiles[c, 1]
                edges_of_type_e = E.chk_profiles[c, 1 + e]
                
                if edges_of_type_e > 0
                    # Sum the incoming variance from ALL edges connected to this check
                    sum_inv_J = 0.0
                    for k in 1:ne
                        edges_of_type_k = E.chk_profiles[c, 1 + k]
                        # Extrinsic principle: pull out ONE edge of the current type 'e'
                        actual_edges = (k == e) ? edges_of_type_k - 1 : edges_of_type_k
                        if actual_edges > 0
                            sum_inv_J += actual_edges * (_inv_J_MI(1.0 - I_V2C[k]))^2
                        end
                    end
                    
                    # Convert sum back to MI and weight it
                    mi_sum += frac * edges_of_type_e * (1.0 - _J_MI(sqrt(sum_inv_J)))
                    weight_sum += frac * edges_of_type_e
                end
            end
            if weight_sum > 0
                new_I_C2V[e] = mi_sum / weight_sum
            end
        end
        I_C2V .= new_I_C2V
        
        # ---------------------------------------------------------
        # 2. Variable Node Update (VND): Combine MI in the variance domain
        # ---------------------------------------------------------
        new_I_V2C = zeros(Float64, ne)
        for e in 1:ne
            mi_sum = 0.0
            weight_sum = 0.0
            
            for v in 1:nv_classes
                frac = E.var_profiles[v, 1]
                is_transmitted = E.var_profiles[v, 2]
                edges_of_type_e = E.var_profiles[v, 2 + e]
                
                if edges_of_type_e > 0
                    # Start with channel variance
                    sum_inv_J = (is_transmitted * sigma_ch)^2
                    
                    # Sum incoming MI from all checks
                    for k in 1:ne
                        edges_of_type_k = E.var_profiles[v, 2 + k]
                        actual_edges = (k == e) ? edges_of_type_k - 1 : edges_of_type_k
                        if actual_edges > 0
                            sum_inv_J += actual_edges * (_inv_J_MI(I_C2V[k]))^2
                        end
                    end
                    
                    mi_sum += frac * edges_of_type_e * _J_MI(sqrt(sum_inv_J))
                    weight_sum += frac * edges_of_type_e
                end
            end
            if weight_sum > 0
                new_I_V2C[e] = mi_sum / weight_sum
            end
        end
        
        # Check convergence
        max_delta = maximum(abs.(new_I_V2C .- I_V2C))
        I_V2C .= new_I_V2C
        
        # If all edge messages reach perfect certainty (1.0)
        if minimum(I_V2C) > 0.999
            return true
        end
        
        # If the solver stalled (hit the error floor or failed the waterfall)
        if max_delta < tol
            return false
        end
    end
    
    return false
end

"""
$(TYPEDSIGNATURES)

Run Multi-Edge Type (MET) Density Evolution using Gaussian Approximation for the BAWGN channel.
"""
function density_evolution_MET_GA(E::METEnsemble, Ch::BAWGNChannel; max_iters::Int=1000, tol::Float64=1e-5)
    # Convert AWGN noise standard deviation to LLR standard deviation
    sigma_ch = 2.0 / Ch.param
    return density_evolution_MET_GA(E, sigma_ch; max_iters=max_iters, tol=tol)
end

"""
$(TYPEDSIGNATURES)

Run Multi-Edge Type (MET) Density Evolution for any arbitrary discrete or continuous channel 
using the AWGN-equivalent capacity approximation.
"""
function density_evolution_MET_GA(E::METEnsemble, Ch::AbstractChannel; max_iters::Int=1000, tol::Float64=1e-5)
    # 1. Compute the exact Shannon capacity of the arbitrary channel
    I_ch = capacity(Ch)
    
    # 2. Map the capacity to an AWGN-equivalent LLR standard deviation
    sigma_ch = _inv_J_MI(I_ch)
    
    # 3. Run the MET GA engine
    return density_evolution_MET_GA(E, sigma_ch; max_iters=max_iters, tol=tol)
end

"""
$(TYPEDSIGNATURES)

Find the optimal AWGN noise threshold (maximum standard deviation `σ_n`) for a 
Multi-Edge Type (MET) ensemble using Gaussian Approximation.
"""
function optimal_threshold(E::METEnsemble, ::Type{BAWGNChannel}; tol::Float64=1e-4)
    low_σ = 0.1   # Low noise -> Should trivially decode
    high_σ = 3.0  # High noise -> Should fail
    
    # Ensure our upper bound actually fails
    while density_evolution_MET_GA(E, BAWGNChannel(high_σ))
        high_σ *= 2.0
    end
    
    # Fast bisection
    for _ in 1:50
        mid_σ = (low_σ + high_σ) / 2.0
        
        if density_evolution_MET_GA(E, BAWGNChannel(mid_σ))
            low_σ = mid_σ  # Decoded! Push the noise higher.
        else
            high_σ = mid_σ # Failed! Lower the noise.
        end
        
        if (high_σ - low_σ) < tol
            break
        end
    end
    
    return low_σ
end

"""
Internal binary search engine. Finds the maximum parameter `param` for which 
`test_func(param)` evaluates to `true`.
"""
function _bisection_threshold(test_func::Function, low::Float64, high::Float64, tol::Float64, expand_high::Bool)
    # If the bounds are open (like AWGN or Fading), dynamically expand the upper bound
    if expand_high
        while test_func(high) && high < 50.0
            high *= 2.0
        end
    end
    
    for _ in 1:100
        mid = (low + high) / 2.0
        if test_func(mid)
            low = mid  # Ensemble decoded! It can handle a worse channel.
        else
            high = mid # Ensemble failed! It needs a better channel.
        end
        
        if (high - low) < tol
            break
        end
    end
    return low
end

# ---------------------------------------------------------
# Evaluation Wrappers
# ---------------------------------------------------------

# For MET Ensembles, we just call the function we already wrote
_decodes_GA(E::METEnsemble, Ch::AbstractChannel) = density_evolution_MET_GA(E, Ch)

# For classical LDPC Ensembles, we map the capacity to the AWGN equivalent σ
function _decodes_GA(E::LDPCEnsemble, Ch::AbstractChannel)
    # 1. Find Shannon capacity of the channel
    I_ch = capacity(Ch)
    # 2. Find AWGN LLR standard deviation
    sigma_llr = _inv_J_MI(I_ch)
    # 3. Map to raw AWGN noise standard deviation (since σ_llr = 2 / σ_noise)
    sigma_noise = 2.0 / sigma_llr
    
    λ_vec = Float64.(coeff.(E.λ, 0:degree(E.λ)))
    ρ_vec = Float64.(coeff.(E.ρ, 0:degree(E.ρ)))
    
    return _density_evolution_GA(λ_vec, ρ_vec, sigma_noise)
end

# ---------------------------------------------------------
# Channel-Specific Dispatches
# ---------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Find the optimal crossover probability threshold (`p`) for the Binary Symmetric Channel.
"""
optimal_threshold(E::AbstractLDPCFamily, ::Type{BinarySymmetricChannel}; tol::Float64=1e-5) = 
    _bisection_threshold(p -> _decodes_GA(E, BinarySymmetricChannel(p)), 0.0, 0.5, tol, false)

"""
$(TYPEDSIGNATURES)

Find the optimal crossover probability threshold (`p`) for the Z-Channel.
"""
optimal_threshold(E::AbstractLDPCFamily, ::Type{ZChannel}; tol::Float64=1e-5) = 
    _bisection_threshold(p -> _decodes_GA(E, ZChannel(p)), 0.0, 1.0, tol, false)

"""
$(TYPEDSIGNATURES)

Find the optimal Ergodic noise threshold (`σ`) for the Rayleigh Fading Channel.
"""
optimal_threshold(E::AbstractLDPCFamily, ::Type{RayleighFadingChannel}; tol::Float64=1e-4) = 
    _bisection_threshold(σ -> _decodes_GA(E, RayleighFadingChannel(σ)), 0.1, 3.0, tol, true)

# (If you wish to overwrite the AWGN threshold we wrote earlier to use this clean engine:)
optimal_threshold(E::AbstractLDPCFamily, ::Type{BAWGNChannel}; tol::Float64=1e-4) = 
    _bisection_threshold(σ -> _decodes_GA(E, BAWGNChannel(σ)), 0.1, 3.0, tol, true)
