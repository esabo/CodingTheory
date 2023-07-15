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

"""
    type(Ch::AbstractClassicalNoiseChannel)

Return the type of the noise channel.
"""
type(Ch::AbstractClassicalNoiseChannel) = typeof(Ch)

#############################
      # setter functions
#############################

#############################
     # general functions
#############################

# TODO: re-evaluate which of these are actually useful
_polyeval(x::Real, vec::Vector{<:Real}) = sum(c * x^(i - 1) for (i, c) in enumerate(vec))
_polyeval(x::Real, f::PolyRingElem) = _polyeval(x, Float64.(coeff.(f, 0:degree(f))))

_dpolyeval(x::Real, vec::Vector{<:Real}) = sum((i - 1) * c * x^(i - 2) for (i, c) in enumerate(vec))
_dpolyeval(x::Real, f::PolyRingElem) = _dpolyeval(x, Float64.(coeff.(f, 0:degree(f))))

_integratepolyeval01(vec::Vector{<:Real}) = sum(c / i for (i, c) in enumerate(vec))
_integratepolyeval01(f::PolyRingElem) = _integratepolyeval01(Float64.(coeff.(f, 0:degree(f))))

_dpoly(vec::Vector{<:Real}) = [vec[i] * (i - 1) for i in 2:length(vec)]
_dpoly(f::PolyRingElem) = derivative(f)

_integratepoly(vec::Vector{T}) where T <: Real = [zero(T); [c / i for (i, c) in enumerate(vec)]]
_integratepoly(f::PolyRingElem) = integral(f)

# TODO: move to utils, check if already there, export
# possibly just go ahead and extend to the nonbinary case then call with a 2 here
_binaryentropy(x::Real) = x * log2(1 / x) + (1 - x) * log2(1 / (1 - x))

# these should all be useful. Note that the QQ version of _computeλρ is necessary. The non-QQ version is for when the poly isn't directly callable (as in RealPoly)
_computeavgdegree(f::PolyRingElem) = inv(sum(coeff(f, i) / (i + 1) for i in 0:degree(f)))
_computeLR(λ::PolyRingElem, ρ::PolyRingElem, lavg, ravg) = (integral(λ) * lavg, integral(ρ) * ravg)
_computeLR(λ::PolyRingElem, ρ::PolyRingElem) = _computeLR(λ, ρ, _computeavgdegree(λ), _computeavgdegree(ρ))
_computeλρ(L::PolyRingElem, R::PolyRingElem) = (derivative(L) / _dpolyeval(1, L), derivative(R) / _dpolyeval(1, R))
_computeλρ(L::QQPolyRingElem, R::QQPolyRingElem) = (derivative(L) / derivative(L)(1), derivative(R) / derivative(R)(1))

# TODO: currently returning L2 squared. Should we sqrt it?
function _L2(p1::Vector{Float64}, p2::Vector{Float64})
    @assert length(p1) == length(p2)
    v = p1 .- p2
    v2 = [sum(v[j] * v[k + 1 - j] for j in max(1, k + 1 - length(v)):min(k, length(v))) for k in 1:2length(v) - 1]
    return _integratepolyeval01(v2)
end

Base.hash(Ch::AbstractClassicalNoiseChannel) = hash(Ch.param, hash(typeof(Ch)))
Base.isequal(Ch1::AbstractClassicalNoiseChannel, Ch2::AbstractClassicalNoiseChannel) = typeof(Ch1) == typeof(Ch2) && Ch1.param == Ch2.param

# function Base.setproperty!(Ch::BAWGNChannel, key, val)
#     if key == :param
#         # TODO: should we throw an error for trying to change this at all and tell the user to create a new channel?
#         setfield!(Ch, :capacity, missing)
#     end
#     setfield!(Ch, key, val)
# end

# TODO: say where this is stored in the tutorial
"""
    densityevolution!(E::LDPCEnsemble, Ch::AbstractClassicalNoiseChannel)

Compute the density evolution of the LDPC ensemble given the noise channel.
"""
function densityevolution!(E::LDPCEnsemble, Ch::AbstractClassicalNoiseChannel)
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
    densityevolution!(E::LDPCEnsemble, Ch::AbstractClassicalNoiseChannel)
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
    Ch ∈ keys(E.densityevo) || densityevolution!(E, Ch)
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

    # @show (1 / α, abs(round(1/α) - 1/α))

    _, x = PolynomialRing(RealField(), :x)
    λ = sum(λi * x^i for (i, λi) in enumerate(λvec))
    ρ = x^round(Int, 1 / α)
    return λ, ρ
end

##############################################
########## still need to edit below ##########
##############################################

# function _findlambdagivenrho(ρvec::Vector{Float64}, ε::Float64, lmax::Int)
#     model = Model(GLPK.Optimizer)
#     @variable(model, λ[1:lmax - 1] >= 0)
#     @constraint(model, sum(λ) == 1)
#     for x in 0:0.01:1
#         @constraint(model, ε * sum(λ[i] * (1 - _polyeval(1 - x, ρvec))^i for i in 1:lmax - 1) - x <= 0)
#     end
#     @constraint(model, ε * _dpolyeval(1, ρvec) * λ[1] <= 1)
#     @objective(model, Max, sum(λ[i] / (i + 1) for i in 1:lmax - 1))
#     optimize!(model)
#     @assert termination_status(model) == MOI.OPTIMAL "Didn't find an optimal point"
#     @assert primal_status(model) == MOI.FEASIBLE_POINT "Didn't optimize to a feasible point"
#     λvec = [0.0; value.(λ)]
#     return λvec
# end

# function _findrhogivenlambda(λvec::Vector{Float64}, ε::Float64, rmax::Int)
#     model = Model(GLPK.Optimizer)
#     @variable(model, ρ[1:rmax - 1] >= 0)
#     @constraint(model, sum(ρ) == 1)
#     for x in 0:0.01:1
#         @constraint(model, 1 - x - sum(ρ[i] * (1 - ε * _polyeval(x, λvec))^i for i in 1:rmax - 1) <= 0)
#     end
#     @constraint(model, ε * sum(i * ρ[i] for i in eachindex(ρ)) * λvec[2] <= 1)
#     @objective(model, Min, sum(ρ[i] / (i + 1) for i in 1:rmax - 1))
#     optimize!(model)
#     @assert termination_status(model) == MOI.OPTIMAL "Didn't find an optimal point"
#     @assert primal_status(model) == MOI.FEASIBLE_POINT "Didn't optimize to a feasible point"
#     ρvec = [0.0; value.(ρ)]
#     return ρvec
# end


# # @assert lmax > 1
# # @assert 0 < ε < 1
# # channel ∈ (:BEC, :BSC, :BIAWGN) || throw(ArgumentError("Channel not yet implemented"))
# # decoder ∈ (:BEC, :A, :SP) || throw(ArgumentError("Decoder not supported"))
# # if channel == :BEC
# #     decoder == :BEC || throw(ArgumentError("The only decoder supported for the BEC channel is :BEC"))
# # end

# # R, x = PolynomialRing(RealField(), :x)

# # type check parameters in these public functions
# function optimallambda(ρ, lmax::Int, realparam::Real, vartype::Symbol)
#     vartype ∈ (:r, :ε) || throw(ArgumentError("Vartype must be :r for target rate or :ε for threshold"))

#     return _optimaldistributions(ρ, :ρ, lmax, realparam, vartype)
# end

# function optimalrho(λ, rmax::Int, realparam::Real, vartype::Symbol)
#     vartype ∈ (:r, :ε) || throw(ArgumentError("Vartype must be :r for target rate or :ε for threshold"))

#     return _optimaldistributions(λ, :λ, rmax, realparam, vartype)
# end

# function _optimaldistributions(poly, polytype::Symbol, varmax::Int, realparam::Real, vartype::Symbol)
#     intpoly = _integratepolyeval01(poly)
#     tolerance = 1e-4
#     if vartype == :r
#         maxiter = 100
#         high = 1.0
#         low = 0.0
#         mid = 0.5
#         if polytype == :ρ
#             sol = _findlambdagivenrho(poly, mid, varmax)
#             solrate = 1 - intpoly / _integratepolyeval01(sol)
#         else
#             sol = _findrhogivenlambda(poly, mid, varmax)
#             solrate = 1 - _integratepolyeval01(sol) / intpoly
#         end
#         Δ = solrate - realparam
#         converged = abs(Δ) <= tolerance
#         count = 0
#         tolerance = 1e-9
#         while abs(Δ) > tolerance && count <= maxiter # make like copy below
#             count += 1
#             Δ > 0 ? (low = mid;) : (high = mid;)
#             mid = (high + low) / 2
#             if polytype == :ρ
#                 sol = _findlambdagivenrho(poly, mid, varmax)
#                 solrate = 1 - intpoly / _integratepolyeval01(sol)
#             else
#                 sol = _findrhogivenlambda(poly, mid, varmax)
#                 solrate = 1 - _integratepolyeval01(sol) / intpoly
#             end
#             Δ = solrate - realparam
#             abs(Δ) > tolerance && (converged = true; break;)
#         end
#         # make into poly?
#         # I think no, because RealPoly doesn't seem like a useful final type to me. Probably just as easy to return a vector and let a different function figure out rationals for a desired n, and that function can return a poly with QQ coeffs. -Ben
#         converged && return sol, solrate, mid
#         error("The LP for $polytype did not converge in $maxiter iterations")
#     else
#         if polytype == :ρ
#             sol = _findlambdagivenrho(poly, realparam, varmax)
#             solrate = 1 - intpoly / _integratepolyeval01(sol)
#         else
#             sol = _findrhogivenlambda(poly, realparam, varmax)
#             solrate = 1 - _integratepolyeval01(sol) / intpoly
#         end
#         # should this also return realparam so that it matches the return above? -Ben
#         return sol, solrate
#     end
# end

# function optimallambdaandrho(lmax::Int, rmax::Int, realparam::Real, vartype::Symbol)
#     vartype ∈ (:r, :ε) || throw(ArgumentError("Vartype must be :r for target rate or :ε for threshold"))

#     tolerance = 1e-9
#     # solve until convergence for each ε, see if rates match, else change ε and repeat
#     if vartype == :r
#         maxiter = 100
#         high = 1.0
#         low = 0.0
#         mid = 0.5

#         # TODO: need some initial guesses here
#         # for now, pick ρ(x) = x^(rmax-1) and find λ from that (i.e., start from regular LDPC assumption). -Ben
#         countinner = 0
#         ρ = zeros(rmax); ρ[end] = 1
#         λ = _findlambdagivenrho(ρ, mid, lmax)
#         ρ = _findrhogivenlambda(λ, mid, rmax)
#         λprev = copy(λ)
#         ρprev = copy(ρ)
#         convergedinner = false
#         while countinner <= maxiter
#             countinner += 1
#             λ = _findlambdagivenrho(ρ, mid, lmax)
#             ρ = _findrhogivenlambda(λ, mid, rmax)
#             normλ = _L2(λ, λprev)
#             normρ = _L2(ρ, ρprev)
#             normλ <= tolerance && normρ <= tolerance && (convergedinner = true; break;)
#             λprev .= λ
#             ρprev .= ρ
#         end
#         if !convergedinner
#             # TODO: better error here, just putting something for now
#             error("inner convergence failed")
#         end
#         solrate = 1 - _integratepolyeval01(ρ) / _integratepolyeval01(ρ)
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
#                 λ = _findlambdagivenrho(ρ, mid, lmax)
#                 ρ = _findrhogivenlambda(λ, mid, rmax)
#                 normλ = _L2(λ, λprev)
#                 normρ = _L2(ρ, ρprev)
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
#         λ = _findlambdagivenrho(ρ, realparam, lmax)
#         ρ = _findrhogivenlambda(λ, realparam, rmax)
#         λprev = copy(λ)
#         ρprev = copy(ρ)
#         convergedinner = false
#         while countinner <= maxiter
#             countinner += 1
#             λ = _findlambdagivenrho(ρ, realparam, lmax)
#             ρ = _findrhogivenlambda(λ, realparam, rmax)
#             normλ = _L2(λ, λprev)
#             normρ = _L2(ρ, ρprev)
#             normλ <= tolerance && normρ <= tolerance && (convergedinner = true; break;)
#             λprev .= λ
#             ρprev .= ρ
#         end
#         convergedinner && (return λ, ρ, 1 - _integratepolyeval01(ρ) / _integratepolyeval01(ρ))
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





