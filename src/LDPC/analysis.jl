# Copyright (c) 2023 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # Classical
#############################

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

# TODO: in tutorial, make sure we explain why we can't pass in L, R here
# TODO: write public conversion functions between all the polynomial types
"""
    LDPCEnsemble(λ::PolyRingElem, ρ::PolyRingElem)

Return the LDPC ensemble determined by the variable degree distribution `λ` and the check
degree distribution `ρ`, both from an edge perspective.
"""
function LDPCEnsemble(λ::PolyRingElem, ρ::PolyRingElem)
    # TODO: check these are proper polynomials, degrees sum to 1, are positive (check if L/R)
    lavg = _computeavgdegree(λ)
    ravg = _computeavgdegree(ρ)
    L, R = _computeLR(λ, ρ, lavg, ravg)
    designrate = Float64(1 - lavg / ravg)
    densityevo = Dict{AbstractClassicalNoiseChannel, NTuple{2, Vector{Float64}}}()
    threshold = Dict{Type, Float64}()
    return LDPCEnsemble(λ, ρ, L, R, Float64(lavg), Float64(ravg), designrate, densityevo, threshold)
end

#############################
      # getter functions
#############################

"""
    erasureprobability(Ch::BinaryErasureChannel)

Return the erasure probability of the binary erasure channel.
"""
erasureprobability(Ch::BinaryErasureChannel) = Ch.param

"""
    crossoverprobability(Ch::BinarySymmetricChannel)

Return the crossover probability of the binary symmetric channel.
"""
crossoverprobability(Ch::BinarySymmetricChannel) = Ch.param

"""
    standarddeviation(Ch::BAWGNChannel)

Return the standard deviation of the BAWGN channel.
"""
standarddeviation(Ch::BAWGNChannel) = Ch.param

"""
    variance(Ch::BAWGNChannel)

Return the variance of the BAWGN channel.
"""
variance(Ch::BAWGNChannel) = Ch.param^2

"""
    capacity(Ch::AbstractClassicalNoiseChannel)

Return the capacity of the noise channel.
"""
capacity(Ch::AbstractClassicalNoiseChannel) = ismissing(Ch.capacity) ? computecapacity(Ch) : Ch.capacity

#############################
      # setter functions
#############################

#############################
     # general functions
#############################

# TODO: are _dpoly(...) and _integratepoly(...) useful? currently unused, probably delete
_dpoly(vec::Vector{<:Real}) = [vec[i] * (i - 1) for i in 2:length(vec)]
_dpoly(f::PolyRingElem) = derivative(f)

_integratepoly(vec::Vector{T}) where T <: Real = [zero(T); [c / i for (i, c) in enumerate(vec)]]
_integratepoly(f::PolyRingElem) = integral(f)

# _polyeval, _dpolyeval, and _integratepoly01 are all used
_polyeval(x::Real, vec::Vector{<:Real}) = sum(c * x^(i - 1) for (i, c) in enumerate(vec))
_polyeval(x::Real, f::PolyRingElem) = _polyeval(x, Float64.(coeff.(f, 0:degree(f))))

_dpolyeval(x::Real, vec::Vector{<:Real}) = sum((i - 1) * c * x^(i - 2) for (i, c) in enumerate(vec))
_dpolyeval(x::Real, f::PolyRingElem) = _dpolyeval(x, Float64.(coeff.(f, 0:degree(f))))

_integratepoly01(vec::Vector{<:Real}) = sum(c / i for (i, c) in enumerate(vec))
_integratepoly01(f::PolyRingElem) = _integratepoly01(Float64.(coeff.(f, 0:degree(f))))

# TODO: move to utils, check if already there, export
# possibly just go ahead and extend to the nonbinary case then call with a 2 here
_binaryentropy(x::Real) = x * log2(1 / x) + (1 - x) * log2(1 / (1 - x))

# these should all be useful. Note that the QQ version of _computeλρ is necessary for the output to also be a QQPoly. The non-QQ version is for when the poly isn't directly callable (as in RealPoly)
_computeavgdegree(f::PolyRingElem) = inv(sum(coeff(f, i) / (i + 1) for i in 0:degree(f)))
_computeLR(λ::PolyRingElem, ρ::PolyRingElem, lavg, ravg) = (integral(λ) * lavg, integral(ρ) * ravg)
_computeLR(λ::PolyRingElem, ρ::PolyRingElem) = _computeLR(λ, ρ, _computeavgdegree(λ), _computeavgdegree(ρ))
_computeλρ(L::PolyRingElem, R::PolyRingElem) = (derivative(L) / _dpolyeval(1, L), derivative(R) / _dpolyeval(1, R))
_computeλρ(L::QQPolyRingElem, R::QQPolyRingElem) = (derivative(L) / derivative(L)(1), derivative(R) / derivative(R)(1))

function _L2distsq(p1::Vector{Float64}, p2::Vector{Float64})
    @assert length(p1) == length(p2)
    v = p1 .- p2
    v2 = [sum(v[j] * v[k + 1 - j] for j in max(1, k + 1 - length(v)):min(k, length(v))) for k in 1:2length(v) - 1]
    return _integratepoly01(v2)
end

Base.hash(Ch::AbstractClassicalNoiseChannel) = hash(Ch.param, hash(typeof(Ch)))
Base.isequal(Ch1::AbstractClassicalNoiseChannel, Ch2::AbstractClassicalNoiseChannel) = typeof(Ch1) == typeof(Ch2) && Ch1.param == Ch2.param

# function Base.setproperty!(Ch::BAWGNChannel, key, val)
#     key == :capacity && (setfield!(Ch, key, val);)
#     key == :capacity || @warn "Channel not updated. Create a new channel instead of changing the noise on an existing channel."
# end

function _densityevolution!(E::LDPCEnsemble, Ch::AbstractClassicalNoiseChannel)
    if isa(Ch, BinaryErasureChannel)
        λvec = Float64.(coeff.(E.λ, 0:degree(E.λ)))
        ρvec = Float64.(coeff.(E.ρ, 0:degree(E.ρ)))
        E.densityevo[Ch] = _densityevolutionBEC(λvec, ρvec, Ch.param)
    else
        error("Only BEC has been implemented so far")
    end
    return nothing
end

"""
    densityevolution(E::LDPCEnsemble, Ch::AbstractClassicalNoiseChannel)

Return the density evolution of the LDPC ensemble given the noise channel.
"""
function densityevolution(E::LDPCEnsemble, Ch::AbstractClassicalNoiseChannel)
    Ch ∈ keys(E.densityevo) || _densityevolution!(E::LDPCEnsemble, Ch::AbstractClassicalNoiseChannel)
    return E.densityevo[Ch]
end

function _densityevolutionBEC(λ::Vector{<:Real}, ρ::Vector{<:Real}, ε::Real;
    maxiters::Int=500, tol::Float64=1e-9)

    iter = 0
    evox = [ε]
    evoy = [1.0]
    while evox[end] > tol && iter < maxiters
        iter += 1
        push!(evoy, 1 - _polyeval(1 - evox[end], ρ))
        push!(evox, ε * _polyeval(evoy[end], λ))
    end
    return evox, evoy
end

"""
    plotEXITchart(E::LDPCEnsemble, Ch::AbstractClassicalNoiseChannel; tol::Float64=1e-9)

Return a plot of the EXIT chart for the ensemble given the channel up to a numerical tolerance of `tol`.
"""
function plotEXITchart(E::LDPCEnsemble, Ch::AbstractClassicalNoiseChannel; tol::Float64=1e-9)
    @assert isa(Ch, BinaryErasureChannel) "Only BEC is implemented so far"
    x = 0:0.01:1
    c = 1 .- [_polyeval(x, E.ρ) for x in 1 .- x]
    v = Ch.param .* [_polyeval(x, E.λ) for x in x]
    Ch ∈ keys(E.densityevo) || _densityevolution!(E, Ch)
    evox, evoy = E.densityevo[Ch]
    ind = findfirst(evox .<= tol)
    title = if isnothing(ind)
        "Failed to converge after $(length(evox) - 1) iterations (tol = $tol)"
    else
        "Converged after $(ind - 1) iterations (tol = $tol), \$\\varepsilon = $(Ch.param)\$"
    end
    p = plot(x, c, label = "\$c(x)\$")
    plot!(p, v, x, label = "\$v_{\\varepsilon}^{-1}(x)\$")
    plot!(p, evox[1:ind], evoy[1:ind], seriestype=:steppre, linestyle=:dash, linecolor=:black, label=false)
    plot!(p, legend=:bottomright, xlims = (0, min(1, Ch.param * 1.2)), ylims = (0, 1), title=title)
    return p
end

# TODO: what else should we accept here and under which do we want to store this and threshold?
"""
    multiplicativegap(E::LDPCEnsemble, Ch::AbstractClassicalNoiseChannel)

Return the multiplicative gap of the ensemble with respect to channel.
"""
# function multiplicativegap(E::LDPCEnsemble, Ch::AbstractClassicalNoiseChannel)
#     threshold = ismissing(___.threshold) ? ___.threshold : computethreshold(E, Ch)
#     return (1 - threshold - E.designrate) / (1 - E.designrate)
# end

"""
    multiplicativegaplowerbound(E::LDPCEnsemble)

Return a lower bound on the multiplicative gap of the ensemble
"""
multiplicativegaplowerbound(E::LDPCEnsemble) = (E.designrate^E.ravg * (1 - E.designrate)) / (1 + E.designrate^E.ravg * (1 - E.designrate))

"""
    densitylowerbound(Ch::AbstractClassicalNoiseChannel, gap::Real)

Return a lower bound on the density of a (full rank) parity-check matrix for the channel
given the multiplicative gap.
"""
function densitylowerbound(Ch::AbstractClassicalNoiseChannel, gap::Real)
    0 < gap < 1 || throw(DomainError("Multiplicative gap should be in (0, 1)"))
    if isa(Ch, BinaryErasureChannel)
        temp = log(1 - Ch.param)
        return (Ch.param * (log(gap) - (log(Ch.param) - temp))) / ((1 - Ch.param) * (1 - gap) * temp)
    else
        @error "Not yet implemented"
    end
end

"""
    checkconcentrateddegreedistribution(Ch::BinaryErasureChannel, gap::Real)

Return the check-concentrated degree distribution `(λ, ρ)` for the binary erasure channel
given the desired multiplicative gap.
"""
function checkconcentrateddegreedistribution(Ch::BinaryErasureChannel, gap::Real)
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

    λvec = zeros(N - 1)
    λvec[1] = α
    for i in 2:N - 1
        λvec[i] = ((i - 1) / i) * (1 - α / (i - 1)) * λvec[i - 1]
    end
    norm = sum(λvec)
    λvec ./= norm

    _, x = PolynomialRing(RealField(), :x)
    λ = sum(λi * x^i for (i, λi) in enumerate(λvec))
    ρ = x^round(Int, 1 / α)
    return λ, ρ
end

function _findlambdagivenrho(ρ::Union{Vector{Float64}, PolyRingElem}, ε::Float64, lmax::Int; Δ = 0.001)
    model = Model(GLPK.Optimizer)
    @variable(model, λ[1:lmax - 1] >= 0)
    @constraint(model, sum(λ) == 1)
    for x in 0:Δ:1
        @constraint(model, ε * sum(λ[i] * (1 - _polyeval(1 - x, ρ))^i for i in 1:lmax - 1) - x <= 0)
    end
    @constraint(model, ε * _dpolyeval(1, ρ) * λ[1] <= 1)
    @objective(model, Max, sum(λ[i] / (i + 1) for i in 1:lmax - 1))
    optimize!(model)
    termination_status(model) == MOI.INFEASIBLE && throw(DomainError("No solution exists"))
    @assert termination_status(model) == MOI.OPTIMAL "Didn't find an optimal point"
    return [0.0; value.(λ)], model
end

function _findrhogivenlambda(λ::Union{Vector{Float64}, PolyRingElem}, ε::Float64, rmax::Int; Δ = 0.001)
    model = Model(GLPK.Optimizer)
    @variable(model, ρ[1:rmax - 1] >= 0)
    @constraint(model, sum(ρ) == 1)
    for x in 0:Δ:1
        @constraint(model, 1 - x - sum(ρ[i] * (1 - ε * _polyeval(x, λ))^i for i in 1:rmax - 1) <= 0)
    end
    temp = isa(λ, PolyRingElem) ? Float64(coeff(λ, 1)) : λ[2]
    @constraint(model, ε * sum(i * ρ[i] for i in eachindex(ρ)) * temp <= 1)
    @objective(model, Min, sum(ρ[i] / (i + 1) for i in 1:rmax - 1))
    optimize!(model)
    termination_status(model) == MOI.INFEASIBLE && throw(DomainError("No solution exists"))
    @assert termination_status(model) == MOI.OPTIMAL "Didn't find an optimal point"
    return [0.0; value.(ρ)], model
end

"""
    optimallambda(ρ, lmax, param, vartype; Δ = 1e-3)

Find the optimal variable node distribution `λ` given the check node
distribution `ρ`, maximum variable node degree `lmax`, and target parameter
`param` which refers to threshold if `vartype == :ε` or rate if `vartype == :r`.

# Notes
* `Δ` refers to the step size for `x` in the LP to solve for `λ`.
"""
function optimallambda(ρ, lmax::Int, param::Float64, vartype::Symbol; Δ = 0.001)
    vartype ∈ (:r, :ε) || throw(ArgumentError("Vartype must be :r for target rate or :ε for threshold"))
    vartype == :r && param >= 1 - 2 * _integratepoly01(ρ) && throw(ArgumentError("This rate is unachieveable with the given ρ."))
    # TODO: check for when vartype == :ε as well

    λvec, r, ε =  _optimaldistributions(ρ, :ρ, lmax, param, vartype, Δλ = 0.001)
    _, x = PolynomialRing(RealField(), :x)
    λ = sum(c * x^(i - 1) for (i, c) in enumerate(λvec))
    return (λ = λ, r = r, ε = ε)
end

"""
    optimalrho(λ, rmax, param, vartype; Δ = 1e-3)

Find the optimal check node distribution `ρ` given the variable node
distribution `λ`, maximum check node degree `rmax`, and target parameter `param`
which refers to threshold if `vartype == :ε` or rate if `vartype == :r`.

# Notes
* `Δ` refers to the step size for x in the LP to solve for ρ.
"""
function optimalrho(λ, rmax::Int, param::Float64, vartype::Symbol, Δ = 1e-3)
    vartype ∈ (:r, :ε) || throw(ArgumentError("Vartype must be :r for target rate or :ε for threshold"))
    vartype == :r && param >= 1 - 1 / (rmax * _integratepoly01(λ)) && throw(ArgumentError("This rate is unachieveable with the given rmax and λ."))
    # TODO: check for when vartype == :ε as well

    ρvec, r, ε = _optimaldistributions(λ, :λ, rmax, param, vartype, Δρ = 0.001)
    _, x = PolynomialRing(RealField(), :x)
    ρ = sum(c * x^(i - 1) for (i, c) in enumerate(ρvec))
    return (ρ = ρ, r = r, ε = ε)
end

function _optimaldistributions(poly, polytype::Symbol, varmax::Int, realparam::Float64, vartype::Symbol; Δλ = 0.001, Δρ = 0.001)
    intpoly = _integratepoly01(poly)
    if vartype == :r
        tolerance = 1e-6
        maxiter = 100
        high = 1.0
        low = 0.0
        mid = 0.5
        count = 0
        Δ = -Inf
        while !(0 <= Δ <= tolerance) && count < maxiter
            count += 1
            mid = (high + low) / 2
            Δ, sol = try
                if polytype == :ρ
                    sol, _ = _findlambdagivenrho(poly, mid, varmax, Δ = Δλ)
                    solrate = 1 - intpoly / _integratepoly01(sol)
                else
                    sol, _ = _findrhogivenlambda(poly, mid, varmax, Δ = Δρ)
                    solrate = 1 - _integratepoly01(sol) / intpoly
                end
                (solrate - realparam, sol)
            catch
                (-Inf, Float64[])
            end
            Δ > 0 ? (low = mid;) : (high = mid;)
        end
        0 <= Δ <= tolerance || error("Solution for $(polytype == :ρ ? :λ : :ρ) did not converge in $maxiter iterations")
        return sol, solrate, mid
    else
        if polytype == :ρ
            sol, _ = _findlambdagivenrho(poly, realparam, varmax, Δ = Δλ)
            solrate = 1 - intpoly / _integratepoly01(sol)
        else
            sol, _ = _findrhogivenlambda(poly, realparam, varmax, Δ = Δρ)
            solrate = 1 - _integratepoly01(sol) / intpoly
        end
        return sol, solrate, realparam
    end
end

function _findlambdaandrho(lmax::Int, rmax::Int, ε::Float64; Δρ = 0.01, Δλ = 0.001)
    # TODO: probably try a smallish Δρ (0.05?) and then bisect in the interval where it's best.
    c = 0:Δρ:1
    iters = length(c)
    ρ = zeros(rmax, iters)
    λ = zeros(lmax, iters)
    rates = fill(-Inf, iters)
    for i in 1:iters
        ρ[end, i] = c[i]
        ρ[end - 1, i] = (1 - c[i])
        try
            λ[:, i] .= _findlambdagivenrho(ρ[:, i], ε, lmax, Δ = Δλ)[1]
            rates[i] = 1 - _integratepoly01(ρ[:, i]) / _integratepoly01(λ[:, i])
        catch end
    end
    i = argmax(rates)
    isfinite(rates[i]) || throw(ArgumentError("No solution for given parameters"))
    return λ[:, i], ρ[:, i], rates[i]
end

"""
    optimallambdaandrho(lmax, rmax, param, vartype; Δρ = 1e-2, Δλ = 1e-3)

Find the optimal distribution pair λ, ρ given the `param`, where `param` is
either a threshold if `vartype == :ε` or a target rate if `vartype == :r`.

# Notes
* `Δρ` gives the step size for possible values of `c` where
    ρ = (1 - c) * x^(rmax - 2) + c * x^(rmax - 1)
* `Δλ` gives the step size for values of `x` in the LP for finding λ given ρ.
"""
function optimallambdaandrho(lmax::Int, rmax::Int, param::Float64, vartype::Symbol; Δρ = 0.01, Δλ = 0.001)
    if vartype == :r
        tolerance = 1e-6
        maxiter = 100
        high = 1.0
        low = 0.0
        mid = 0.5
        count = 0
        Δ = -Inf
        while !(0 <= Δ <= tolerance) && count < maxiter
            count += 1
            mid = (high + low) / 2
            Δ, λvec, ρvec = try
                λ, ρ, solrate = _findlambdaandrho(lmax, rmax, mid, Δρ = Δρ, Δλ = Δλ)
                (solrate - param, λ, ρ)
            catch
                (-Inf, Float64[], Float64[])
            end
            Δ > 0 ? (low = mid;) : (high = mid;)
        end
        0 <= Δ <= tolerance || error("Solution for $(polytype == :ρ ? :λ : :ρ) did not converge in $maxiter iterations")
        _, x = PolynomialRing(RealField(), :x)
        λ = sum(c * x^(i - 1) for (i, c) in enumerate(λvec))
        ρ = sum(c * x^(i - 1) for (i, c) in enumerate(ρvec))
        return (λ = λ, ρ = ρ, r = solrate, ε = mid)
    elseif vartype == :ε
        λvec, ρvec, solrate = _findlambdaandrho(lmax, rmax, param, Δρ = Δρ, Δλ = Δλ)
        _, x = PolynomialRing(RealField(), :x)
        λ = sum(c * x^(i - 1) for (i, c) in enumerate(λvec))
        ρ = sum(c * x^(i - 1) for (i, c) in enumerate(ρvec))
        return (λ = λ, ρ = ρ, r = solrate, ε = param)
    end
    throw(ArgumentError("vartype must be :r or :ε"))
end

"""
    optimalthreshold(λ, ρ; Δ = 1e-4)

Given distributions λ and ρ, find the optimal threshold under BP.

# Notes
* For checking stability, it can be useful to use, e.g., `Δ = BigFloat("1e-7")`
"""
function optimalthreshold(λ::Union{Vector{<:Real}, PolyRingElem}, ρ::Union{Vector{<:Real}, PolyRingElem}; Δ::T = 1e-4) where T <: Real
    xs = Δ:Δ:one(T)
    minimum(x / (_polyeval(1 - _polyeval(1 - x, ρ), λ)) for x in xs)
end

convolution(v::Vector{<:Real}...) = real.(ifft(reduce(.*, map(fft, v))))
convolution(v::Vector{<:Real}, num::Int) = real.(ifft(fft(v) .^ num))

function testconv(v::Vector{T}, x::Vector{T}, w::Vector{T} = ones(T, length(v))) where T <: Real
    @assert length(v) == length(w) == length(x)
    Δx = x[2] - x[1]
    @assert all(x[i] - x[i-1] ≈ Δx for i in 2:length(x))
    k = round(Int, 1 - x[1] / Δx)
    temp = w .* v
    [sum(temp[j] * v[k+i-j] for j in eachindex(v) if 1 <= k + i - j <= length(v)) for i in eachindex(v)]
end

function testconv2(v::Vector{T}, x::Vector{T}, w::Vector{T} = ones(T, length(v))) where T <: Real
    @assert length(v) == length(w) == length(x)
    Δx = x[2] - x[1]
    @assert all(x[i] - x[i-1] ≈ Δx for i in 2:length(x))
    k = round(Int, 1 - x[1] / Δx)
    paddedv = vcat(zeros(T, length(v)), v, zeros(T, length(v)))
    paddedconv = real.(ifft(fft(paddedv) .^ 2))
    circshift(paddedconv, k + length(x))[length(v)+1:2length(v)] .* w
end

function _LdensitytoDdensity(a::Vector{<:Real}, x::Vector{<:Real})
    @assert length(a) == length(x)
    xd = tanh.(x ./ 2)
end

function _DdensitytoGdensity(a::Vector{<:Real}, x::Vector{<:Real})
    @assert length(a) == length(x)
end

function _GdensitytoDdensity(a::Vector{<:Real}, x::Vector{<:Real})
    @assert length(a) == length(x)
end

function _DdensitytoLdensity(a::Vector{<:Real}, x::Vector{<:Real})
    @assert length(a) == length(x)
end

function _LdensitytoGdensity(a::Vector{<:Real}, x::Vector{<:Real})
    @assert length(a) == length(x)
end

function _GdensitytoLdensity(a::Vector{<:Real}, x::Vector{<:Real})
    @assert length(a) == length(x)
end

# This one isn't used for density evolution, might be useful later though
function Gconvolution(x::Vector{<:Real}, v::Vector{T}...) where T <: Real
    vg = Vector{T}[]
    temp , xg = _LdensitytoGdensity(v[1], x)
    push!(vg, temp)
    for i in 2:length(v)
        push!(vg, _LdensitytoGdensity(v[i], x)[1])
    end
    return _GdensitytoLdensity(vg, xg)[1]
end

# this one is used for density evolution
function Gconvolution(x::Vector{<:Real}, v::Vector{<:Real}, num::Int)
    vg, xg = _LdensitytoGdensity(v, x)
    return _GdensitytoLdensity(convolution(vg, num), xg)[1]
end

# polynomial evaluation of a density
function _polyeval(a::Vector{<:Real}, vec::Vector{<:Real}, op::Function)
    sum(c * op(a, i - 1) for (i, c) in enumerate(vec))
end
function _polyeval(a::Vector{<:Real}, f::PolyRingElem, op::Function)
    _polyeval(a, Float64.(coeff.(f, 0:degree(f))), op)
end

function _densityevolutionBMS(λ::Vector{<:Real}, ρ::Vector{<:Real}, a::Vector{T},
                              x::Vector{<:Real}; maxiters::Int=10) where T <: Real

    @assert length(a) == length(x)
    # also assert that x is evenly spaced?

    iter = 0
    evoa = Vector{T}[]
    evob = Vector{T}[]
    while iter < maxiters
        iter += 1
        push!(evob, _polyeval(iter == 1 ? a : evoa[end], ρ, (a, i) -> Gconvolution(x, a, i)))
        push!(evoa, _polyeval(evob[end], λ, convolution))
    end
    return evoa, evob
end

struct _LDensity
    N::Int
    δ::Float64
    data::Vector{Float64}
    # `data` is a vector of length 2N+1 on the quadrature -δN:δ:δN. It represents a probability dist:
    # 0 <= sum(data) * δ <= 1, where the remaining density is assumed to be at ∞
end
_LDensity(N::Int, δ::Float64) = _LDensity(N, δ, zeros(2N+1))

function _densityevolutionBMS(λ::Vector{<:Real}, ρ::Vector{<:Real}, initial::_LDensity; maxiters::Int=10)
    N = initial.N
    δ = initial.δ

    # test values, delete later
    N = 40
    δ = 0.05

    # used to calculate b
    Q(i::Int, j::Int) = round(Int, ((i*δ) ⊞ (j*δ)) / δ)
    TQ(i::Int, k::Int) = k == -1 ? N + 1 : minimum(j for j in i:10N if Q(i, j) >= i - k; init = N + 1)
    TQ_precomputed = [TQ(i, k) for i in 0:N, k in -1:ceil(Int, log(2)/δ - 0.5)]

    t = ceil(Int, log2(3N+3))
    afft = zeros(ComplexF64, 2^t)
    temp = exp.(.-range(0, δ * N, length = N + 1) ./ 2)
    initialfft = copy(afft)
    initialfft[1:N + 1] .= initial.data[N + 1:end] .* temp
    fft!(initialfft)

    iter = 0
    a = deepcopy(initial)
    evoa = _LDensity[]
    evob = _LDensity[]
    while iter < maxiters
        # `initial` is the log-likelihood distribution of the channel message
        # `b` is the log-likelihood distribution after CN processing
        # `a` is the log-likelihood distribution after VN processing
        iter += 1

        b = _LDensity(N, δ)
        bp = a.data[N + 1:end] .+ a.data[N + 1:-1:1]
        bp[1] = a.data[N + 1]
        bm = a.data[N + 1:end] .- a.data[N + 1:-1:1]
        Bp = (1 - sum(a.data) * δ) .+ [sum(bp[j] for j in i:length(bp)) for i in eachindex(bp)]
        push!(Bp, 0.0)
        Bm = (1 - sum(a.data) * δ) .+ [sum(bm[j] for j in i:length(bm)) for i in eachindex(bm)]
        push!(Bm, 0.0)
        bptemp = copy(bp)
        bmtemp = copy(bm)
        btotal = zeros(ComplexF64, 2^t)
        for i in 2:length(ρ)
            ;
        end
        push!(evob, b)

        a = _LDensity(N, δ)
        afft[1:N + 1] .= b.data[N + 1:end] .* temp
        afft[N + 2:end] .= zero(ComplexF64)
        fft!(afft)
        atemp = copy(ac)
        atotal = zeros(ComplexF64, 2^t)
        for i in 2:length(λ)
            atemp .*= ac
            atotal .+= λ[i] .* atemp
        end
        atotal .*= initialfft
        ifft!(atotal)
        a.data[N + 1:end] .+= real.(atotal[1:N + 1]) .* exp.(0:δ:N*δ ./ 2) ./ δ
        a.data[1:N] .= a.data[end:-1:N + 2] .* exp.(-N*δ:δ:-δ)
        push!(evoa, a)
    end
    return evoa, evob
end


#########################################################################
########## below is old code, delete after verifying the above ##########
#########################################################################

# # @assert lmax > 1
# # @assert 0 < ε < 1
# # channel ∈ (:BEC, :BSC, :BIAWGN) || throw(ArgumentError("Channel not yet implemented"))
# # decoder ∈ (:BEC, :A, :SP) || throw(ArgumentError("Decoder not supported"))
# # if channel == :BEC
# #     decoder == :BEC || throw(ArgumentError("The only decoder supported for the BEC channel is :BEC"))
# # end

# # R, x = PolynomialRing(RealField(), :x)

# function optimallambdaandrho(lmax::Int, rmax::Int, realparam::Float64, vartype::Symbol)
#     vartype ∈ (:r, :ε) || throw(ArgumentError("Vartype must be :r for target rate or :ε for threshold"))
#     vartype == :r && realparam >= 1 - 2/rmax && throw(ArgumentError("This rate is unachieveable with the given rmax."))
#     # TODO: check for when vartype == :ε as well

#     tolerance = 1e-9

#     # initial guess: ρ(x) = x^(rmax - 1)
#     ρ = zeros(rmax); ρ[end] = 1
#     # this makes sense to me, but we could choose for some c ∈ (0, 1]
#     #   ρ(x) = (1 - c) * x^(rmax - 2) + c * x^(rmax - 1)
#     # which would be given by the code:
#     # ρ = zeros(rmax); c = 0.5; ρ[end - 1] = 1 - c; ρ[end] = c;

#     # if we need an initial λ, this would be it:
#     # λ, _, _ = _optimaldistributions(ρ, :ρ, lmax, realparam, vartype)

#     # solve until convergence for each ε, see if rates match, else change ε and repeat
#     if vartype == :r
#         maxiter = 100
#         high = 1.0
#         low = 0.0
#         mid = 0.5

#         countinner = 0
#         λ, _ = _findlambdagivenrho(ρ, mid, lmax)
#         ρ, _ = _findrhogivenlambda(λ, mid, rmax)
#         λprev = copy(λ)
#         ρprev = copy(ρ)
#         convergedinner = false
#         while countinner <= maxiter
#             countinner += 1
#             λ, _ = _findlambdagivenrho(ρ, mid, lmax)
#             ρ, _ = _findrhogivenlambda(λ, mid, rmax)
#             normλ = _L2distsq(λ, λprev)
#             normρ = _L2distsq(ρ, ρprev)
#             normλ <= tolerance && normρ <= tolerance && (convergedinner = true; break;)
#             λprev .= λ
#             ρprev .= ρ
#         end
#         if !convergedinner
#             # TODO: better error here, just putting something for now
#             error("inner convergence failed")
#         end
#         solrate = 1 - _integratepoly01(ρ) / _integratepoly01(λ)
#         Δ = solrate - realparam

#         converged = abs(Δ) <= tolerance
#         count = 0
#         while count <= maxiter && !converged
#             count += 1
#             Δ > 0 ? (low = mid;) : (high = mid;)
#             mid = (high + low) / 2

#             countinner = 0
#             convergedinner = false
#             while countinner <= maxiter
#                 countinner += 1
#                 λ, _ = _findlambdagivenrho(ρ, mid, lmax)
#                 ρ, _ = _findrhogivenlambda(λ, mid, rmax)
#                 normλ = _L2distsq(λ, λprev)
#                 normρ = _L2distsq(ρ, ρprev)
#                 normλ <= tolerance && normρ <= tolerance && (convergedinner = true; break;)
#                 λprev .= λ
#                 ρprev .= ρ
#             end
#             if !convergedinner
#                 # TODO: better error here, just putting something for now
#                 error("inner convergence failed")
#             end
#             converged = abs(Δ) <= tolerance
#         end
#         # TODO: better error here, just putting something for now
#         converged ? (return λ, ρ, mid;) : error("outer convergence failed")
#     else # vartype == :ε
#         countinner = 0
#         λ, _ = _findlambdagivenrho(ρ, realparam, lmax)
#         ρ, _ = _findrhogivenlambda(λ, realparam, rmax)
#         λprev = copy(λ)
#         ρprev = copy(ρ)
#         convergedinner = false
#         while countinner <= maxiter
#             countinner += 1
#             λ, _ = _findlambdagivenrho(ρ, realparam, lmax)
#             ρ, _ = _findrhogivenlambda(λ, realparam, rmax)
#             normλ = _L2distsq(λ, λprev)
#             normρ = _L2distsq(ρ, ρprev)
#             normλ <= tolerance && normρ <= tolerance && (convergedinner = true; break;)
#             λprev .= λ
#             ρprev .= ρ
#         end
#         convergedinner && (return λ, ρ, 1 - _integratepoly01(ρ) / _integratepoly01(λ))
#     end
# end

# function optimalthreshold(λ, ρ)
#     high = 1.0
#     low = 0.0
#     Δ = 1
#     tolerance = 1e-9
#     while Δ > tolerance
#         mid = (high + low) / 2
#         # some way to evaluate f(x) on this range
#         # for x in 0.001:001:1
#         #     # lmax
#         #     mid * sum(λ[i] * (1 - _polyeval(1 - x, ρvec))^i for i in 1:lmax - 1) - x
#         # end

#         # some loop
#         f(low) * f(high) < 0 ? (high = mid;) : (low = mid;)
#         Δ = abs(f(mid))
#     end
#     return mid
# end

#############################
         # Quantum
#############################

#############################
        # constructors
#############################

#############################
      # getter functions
#############################

#############################
      # setter functions
#############################

#############################
     # general functions
#############################





