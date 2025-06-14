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
    density_evo = Dict{AbstractClassicalNoiseChannel,NTuple{2,Vector{Float64}}}()
    threshold = Dict{Type,Float64}()
    return LDPCEnsemble(
        λ,
        ρ,
        L,
        R,
        Float64(l_avg),
        Float64(r_avg),
        design_rate,
        density_evo,
        threshold,
    )
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
_d_poly(vec::Vector{<:Real}) = [vec[i] * (i - 1) for i = 2:length(vec)]
_d_poly(f::PolyRingElem) = derivative(f)

_integrate_poly(vec::Vector{T}) where {T<:Real} =
    [zero(T); [c / i for (i, c) in enumerate(vec)]]
_integrate_poly(f::PolyRingElem) = integral(f)

# _poly_eval, _d_poly_eval, and _integrate_poly_0_1 are all used
_poly_eval(x::Real, vec::Vector{<:Real}) = sum(c * x^(i - 1) for (i, c) in enumerate(vec))
_poly_eval(x::Real, f::PolyRingElem) = _poly_eval(x, Float64.(coeff.(f, 0:degree(f))))

_d_poly_eval(x::Real, vec::Vector{<:Real}) =
    sum((i - 1) * c * x^(i - 2) for (i, c) in enumerate(vec))
_d_poly_eval(x::Real, f::PolyRingElem) = _d_poly_eval(x, Float64.(coeff.(f, 0:degree(f))))

_integrate_poly_0_1(vec::Vector{<:Real}) = sum(c / i for (i, c) in enumerate(vec))
_integrate_poly_0_1(f::PolyRingElem) = _integrate_poly_0_1(Float64.(coeff.(f, 0:degree(f))))

# TODO: move to utils, check if already there, export
# possibly just go ahead and extend to the nonbinary case then call with a 2 here
_binary_entropy(x::Real) = x * log2(1 / x) + (1 - x) * log2(1 / (1 - x))

# these should all be useful. Note that the QQ version of _compute_λ_ρ is necessary for the output to also be a QQPoly. The non-QQ version is for when the poly isn't directly callable (as in RealPoly)
_compute_avg_degree(f::PolyRingElem) = inv(sum(coeff(f, i) / (i + 1) for i = 0:degree(f)))
_compute_L_R(λ::PolyRingElem, ρ::PolyRingElem, l_avg, r_avg) =
    (integral(λ) * l_avg, integral(ρ) * r_avg)
_compute_L_R(λ::PolyRingElem, ρ::PolyRingElem) =
    _compute_L_R(λ, ρ, _compute_avg_degree(λ), _compute_avg_degree(ρ))
_compute_λ_ρ(L::PolyRingElem, R::PolyRingElem) =
    (derivative(L) / _d_poly_eval(1, L), derivative(R) / _d_poly_eval(1, R))
_compute_λ_ρ(L::QQPolyRingElem, R::QQPolyRingElem) =
    (derivative(L) / derivative(L)(1), derivative(R) / derivative(R)(1))

function _L2_dist_sq(p1::Vector{Float64}, p2::Vector{Float64})
    @assert length(p1) == length(p2)
    v = p1 .- p2
    v2 = [
        sum(v[j] * v[k+1-j] for j = max(1, k+1-length(v)):min(k, length(v))) for
        k = 1:(2length(v)-1)
    ]
    return _integrate_poly_0_1(v2)
end

Base.hash(Ch::AbstractClassicalNoiseChannel) = hash(Ch.param, hash(typeof(Ch)))
Base.isequal(Ch1::AbstractClassicalNoiseChannel, Ch2::AbstractClassicalNoiseChannel) =
    typeof(Ch1) == typeof(Ch2) && Ch1.param == Ch2.param

# function Base.setproperty!(Ch::BAWGNChannel, key, val)
#     key == :capacity && (setfield!(Ch, key, val);)
#     key == :capacity || @warn "Channel not updated. Create a new channel instead of changing the noise on an existing channel."
# end

function _density_evolution!(E::LDPCEnsemble, Ch::AbstractClassicalNoiseChannel)
    if isa(Ch, BinaryErasureChannel)
        λ_vec = Float64.(coeff.(E.λ, 0:degree(E.λ)))
        ρ_vec = Float64.(coeff.(E.ρ, 0:degree(E.ρ)))
        E.density_evo[Ch] = _density_evolution_BEC(λ_vec, ρ_vec, Ch.param)
    else
        error("Only BEC has been implemented so far")
    end
    return nothing
end

"""
    density_evolution(E::LDPCEnsemble, Ch::AbstractClassicalNoiseChannel)

Return the density evolution of the LDPC ensemble given the noise channel.
"""
function density_evolution(E::LDPCEnsemble, Ch::AbstractClassicalNoiseChannel)
    Ch ∈ keys(E.density_evo) ||
        _density_evolution!(E::LDPCEnsemble, Ch::AbstractClassicalNoiseChannel)
    return E.density_evo[Ch]
end

function _density_evolution_BEC(
    λ::Vector{<:Real},
    ρ::Vector{<:Real},
    ε::Real;
    max_iters::Int = 500,
    tol::Float64 = 1e-9,
)

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
    EXIT_chart_plot(E::LDPCEnsemble, Ch::AbstractClassicalNoiseChannel; tol::Float64 = 1e-9)

Return a plot of the EXIT chart for the ensemble given the channel up to a numerical tolerance of `tol`.

# Note
- Run `using Makie` to activate this extension.
"""
function EXIT_chart_plot end

# TODO: what else should we accept here and under which do we want to store this and threshold?
"""
    multiplicative_gap(E::LDPCEnsemble, Ch::AbstractClassicalNoiseChannel)

Return the multiplicative gap of the ensemble with respect to channel.
"""
# function multiplicative_gap(E::LDPCEnsemble, Ch::AbstractClassicalNoiseChannel)
#     threshold = ismissing(___.threshold) ? ___.threshold : computethreshold(E, Ch)
#     return (1 - threshold - E.design_rate) / (1 - E.design_rate)
# end

"""
    multiplicative_gap_lower_bound(E::LDPCEnsemble)

Return a lower bound on the multiplicative gap of the ensemble
"""
multiplicative_gap_lower_bound(E::LDPCEnsemble) =
    (E.design_rate^E.r_avg * (1 - E.design_rate)) /
    (1 + E.design_rate^E.r_avg * (1 - E.design_rate))

"""
    density_lower_bound(Ch::AbstractClassicalNoiseChannel, gap::Real)

Return a lower bound on the density of a (full rank) parity-check matrix for the channel
given the multiplicative gap.
"""
function density_lower_bound(Ch::AbstractClassicalNoiseChannel, gap::Real)
    0 < gap < 1 || throw(DomainError("Multiplicative gap should be in (0, 1)"))
    if isa(Ch, BinaryErasureChannel)
        temp = log(1 - Ch.param)
        return (Ch.param * (log(gap) - (log(Ch.param) - temp))) /
               ((1 - Ch.param) * (1 - gap) * temp)
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
    for i = 2:(N-1)
        λ_vec[i] = ((i - 1) / i) * (1 - α / (i - 1)) * λ_vec[i-1]
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


"""
    optimal_threshold(λ, ρ; Δ = 1e-4)

Given distributions λ and ρ, find the optimal threshold under BP.

# Notes
* For checking stability, it can be useful to use, e.g., `Δ = BigFloat("1e-7")`
"""
function optimal_threshold(
    λ::Union{Vector{<:Real},PolyRingElem},
    ρ::Union{Vector{<:Real},PolyRingElem};
    Δ::T = 1e-4,
) where {T<:Real}

    xs = Δ:Δ:one(T)
    minimum(x / (_poly_eval(1 - _poly_eval(1 - x, ρ), λ)) for x in xs)
end

# convolution(v::Vector{<:Real}...) = real.(ifft(reduce(.*, map(fft, v))))
# convolution(v::Vector{<:Real}, num::Int) = real.(ifft(fft(v) .^ num))

# function testconv(v::Vector{T}, x::Vector{T}, w::Vector{T} = ones(T, length(v))) where T <: Real
#     @assert length(v) == length(w) == length(x)
#     Δx = x[2] - x[1]
#     @assert all(x[i] - x[i-1] ≈ Δx for i in 2:length(x))
#     k = round(Int, 1 - x[1] / Δx)
#     temp = w .* v
#     [sum(temp[j] * v[k+i-j] for j in eachindex(v) if 1 <= k + i - j <= length(v)) for i in eachindex(v)]
# end

# function testconv2(v::Vector{T}, x::Vector{T}, w::Vector{T} = ones(T, length(v))) where T <: Real
#     @assert length(v) == length(w) == length(x)
#     Δx = x[2] - x[1]
#     @assert all(x[i] - x[i-1] ≈ Δx for i in 2:length(x))
#     k = round(Int, 1 - x[1] / Δx)
#     paddedv = vcat(zeros(T, length(v)), v, zeros(T, length(v)))
#     paddedconv = real.(ifft(fft(paddedv) .^ 2))
#     circshift(paddedconv, k + length(x))[length(v)+1:2length(v)] .* w
# end

# using FFTW

# struct _L_Density
#     # `data` is a vector of length 2N+1 on the quadrature -δN:δ:δN. It represents a probability dist:
#     # 0 <= sum(data) * δ <= 1, where the remaining density is assumed to be at ∞
#     N::Int
#     δ::Float64
#     data::Vector{Float64}
# end

# _L_Density(N::Int, δ::Float64) = _L_Density(N, δ, zeros(2N+1))

# _DE_quantizer(i::Int, j::Int, δ::Real) = round(Int, ((i * δ) ⊞ (j * δ)) / δ)

# function _density_evolution_BMS(λ::Vector{<:Real}, ρ::Vector{<:Real}, initial::_L_Density; max_iters::Int=10)
#     # set up `initial_fft` just once so it can be reused throughout
#     N = initial.N
#     δ = initial.δ
#     t = ceil(Int, log2(3N + 3)) + 1
#     initial_fft = zeros(ComplexF64, 2^t)
#     initial_fft[1:N + 1] .= initial.data[N + 1:end] .* exp.(.-collect(0:N) .* δ ./ 2)
#     initial_fft[end - N + 1:end] .= initial.data[1:N] .* exp.(.-collect(-N:-1) .* δ ./ 2)
#     fft!(initial_fft)

#     # main loop
#     iter = 0
#     a = initial
#     evo_a = CodingTheory._L_Density[]
#     evo_b = CodingTheory._L_Density[]
#     while iter < max_iters
#         iter += 1

#         b = _check_node_update_BMS_DE_quantized(ρ, a)
#         # b = _check_node_update_BMS_DE(ρ, a)

#         # if the total probability exceeds 1, normalize (this shouldn't be necessary but is mathematically fine)
#         sum(b.data) * δ > 1 && (b.data ./= sum(b.data) * δ;)

#         # a = _variable_node_update_BMS_DE(λ, b, initial_fft)
#         a = _variable_node_update_BMS_DE_bad(λ, b, initial)

#         # if the total probability exceeds 1, normalize (this shouldn't be necessary but is mathematically fine)
#         sum(a.data) * δ > 1 && (a.data ./= sum(a.data) * δ;)

#         # normalize to 1 probability no matter what (this is mathematically incorrect)
#         # a.data ./= sum(a.data) * δ

#         push!(evo_b, b)
#         push!(evo_a, a)

#         # just to easily track how long things are taking while debugging
#         iter % 10 == 0 ? print(iter) : print(".")
#     end
#     print("\n")
#     return evo_a, evo_b
# end

# function _check_node_update_BMS_DE_quantized(ρ::Vector{<:Real}, a::_L_Density)
#     N = a.N
#     δ = a.δ
#     ap = a.data[N + 1:end] .+ a.data[N + 1:-1:1]
#     ap[1] = a.data[N + 1]
#     am = a.data[N + 1:end] .- a.data[N + 1:-1:1]
#     bp = copy(ap)
#     bm = copy(am)
#     total_p = zeros(N + 1)
#     total_m = zeros(N + 1)
#     ainf = 1 - sum(ap) * δ
#     binf = 1 - sum(bp) * δ
#     for l in 2:length(ρ)
#         # update polynomial evaluation
#         total_p .+= bp * ρ[l]
#         total_m .+= bm * ρ[l]

#         # update convolution
#         cp = zeros(N + 1)
#         cm = zeros(N + 1)

#         for i in 1:N + 1
#             k = _DE_quantizer(i - 1, i - 1, δ) + 1
#             cp[k] += ap[i] * bp[i]
#             cm[k] += am[i] * bm[i]
#             for j in i + 1:N + 1
#                 k = _DE_quantizer(i - 1, j - 1, δ) + 1
#                 cp[k] += ap[i] * bp[j] + ap[j] * bp[i]
#                 cm[k] += am[i] * bm[j] + am[j] * bm[i]
#             end
#         end
#         @. cp += ap * binf + ainf * bp
#         @. cm += am * binf + ainf * bm

#         # update bp, bm, binf
#         @. bp = cp * δ
#         @. bm = cm * δ
#         binf = 1 - sum(bp) * δ
#     end

#     b = _L_Density(N, δ)
#     b.data[N + 1] = total_p[1]
#     b.data[N + 2:end] .= (total_p[2:end] .+ total_m[2:end]) ./ 2
#     b.data[1:N] .= (total_p[end:-1:2] .- total_m[end:-1:2]) ./ 2

#     return b
# end

# function _check_node_update_BMS_DE(ρ::Vector{<:Real}, a::_L_Density)
# end

# function _variable_node_update_BMS_DE(λ::Vector{<:Real}, dist::_L_Density, initial_fft::Vector{<:Complex})
#     N = dist.N
#     δ = dist.δ
#     a_fft = zeros(ComplexF64, length(initial_fft))
#     a_fft[1:N + 1] .= dist.data[N + 1:end] .* exp.(.-collect(0:N) .* δ ./ 2)
#     a_fft[end - N + 1:end] .= dist.data[1:N] .* exp.(.-collect(-N:-1) .* δ ./ 2)
#     fft!(a_fft)

#     # need to convolve with itself multiple times, so save the FFT in a_temp
#     a_temp = copy(a_fft)

#     # keep the running total in a_total
#     a_total = zeros(ComplexF64, length(initial_fft))

#     # collect the FFT domain result (still transformed to take advantage of L-symmetry) of the
#     # polynomial λ applied to `dist`
#     for i in 2:length(λ)
#         @. a_total += λ[i] * a_temp
#         @. a_temp *= a_fft * δ
#     end

#     # convolution with the original channel message
#     a_total .*= initial_fft * δ

#     # get back to the proper domain and undo the transformation from above
#     ifft!(a_total)
#     result = _L_Density(N, δ)
#     result.data[N + 1:end] .= real.(a_total[1:N + 1]) .* exp.(collect(0:N) .* δ ./ 2)
#     result.data[1:N] .= real.(a_total[end - N + 1:end]) .* exp.(collect(-N:-1) .* δ ./ 2)

#     return result
# end

# function _variable_node_update_BMS_DE_bad(λ::Vector{<:Real}, b::_L_Density, initial::_L_Density)
#     _conv(x, y) = b.δ * [sum(x[i - k + b.N + 1] * y[k] for k in eachindex(y) if 1 <= i - k + b.N + 1 <= length(x)) for i in eachindex(y)]
#     temp = copy(b.data)
#     a = _L_Density(b.N, b.δ)
#     for i in 2:length(λ)
#         a.data .+= λ[i] * temp
#         temp .= _conv(temp, b.data)
#     end
#     a.data .= _conv(a.data, initial.data)
#     return a
# end

# function DEBMStest()
#     δ = 0.01
#     N = round(Int, 50 / δ)
#     initial = CodingTheory._L_Density(N, δ)

#     # example 4.100, p221
#     sigma = 0.93
#     initial.data .= [(sigma / √(8π)) * exp(-(y - (2 / sigma^2))^2 * sigma^2 / 8) for y in -δ*N:δ:δ*N]
#     λ = [0, 0.212332, 0.197596, 0, 0.0142733, 0.0744898, 0.0379457, 0.0693008, 0.086264, 0, 0.00788586, 0.0168657, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.283047]
#     ρ = [0, 0, 0, 0, 0, 0, 0, 0, 1.0]

#     inds = (1, 5, 10, 14, 15)
#     # inds = (1, 5, 10, 25, 50)
#     # inds = (1, 5, 10, 50, 140)

#     evo_a, evo_b = _density_evolution_BMS(λ, ρ, initial; max_iters = maximum(inds) + 1)

#     plot_delta = 0.5
#     skip = round(Int, plot_delta / δ)
#     x = -δ * N:plot_delta:δ * N
#     y_inds = 1:skip:length(-plot_delta * N:plot_delta:plot_delta * N)

#     # This plot should look like fig 4.101, top left panel
#     plt1 = plot(x, initial.data[y_inds],
#                 title = "\$a_0\$",
#                 xlims = (-10, 45),
#                 ylims = (0,0.25),
#                 legend = false,
#                 framestyle = :box,
#                 yticks = ([0.05, 0.1, 0.15, 0.2], ["0.05", "0.10", "0.15", "0.20"]),
#                 xticks = (-10:5:40, ["-10", "", "0", "", "10", "", "20", "", "30", "", "40"])
#                 )

#     # This should look like fig 4.101, top right panel
#     plt2 = plot(x, evo_b[1].data[y_inds],
#                 title = "\$b_1\$",
#                 xlims = (-10, 45),
#                 ylims = (0,0.25),
#                 legend = false,
#                 framestyle = :box,
#                 yticks = ([0.05, 0.1, 0.15, 0.2], ["0.05", "0.10", "0.15", "0.20"]),
#                 xticks = (-10:5:40, ["-10", "", "0", "", "10", "", "20", "", "30", "", "40"])
#                 )

#     plt3 = plot(x, evo_a[inds[2]].data[y_inds],
#                 title = "\$a_{$(inds[2])}\$",
#                 xlims = (-10, 45),
#                 ylims = (0,0.25),
#                 legend = false,
#                 framestyle = :box,
#                 yticks = ([0.05, 0.1, 0.15, 0.2], ["0.05", "0.10", "0.15", "0.20"]),
#                 xticks = (-10:5:40, ["-10", "", "0", "", "10", "", "20", "", "30", "", "40"])
#                 )
#     plt4 = plot(x, evo_b[inds[2]+1].data[y_inds],
#                 title = "\$b_{$(inds[2]+1)}\$",
#                 xlims = (-10, 45),
#                 ylims = (0,0.25),
#                 legend = false,
#                 framestyle = :box,
#                 yticks = ([0.05, 0.1, 0.15, 0.2], ["0.05", "0.10", "0.15", "0.20"]),
#                 xticks = (-10:5:40, ["-10", "", "0", "", "10", "", "20", "", "30", "", "40"])
#                 )
#     plt5 = plot(x, evo_a[inds[3]].data[y_inds],
#                 title = "\$a_{$(inds[3])}\$",
#                 xlims = (-10, 45),
#                 ylims = (0,0.25),
#                 legend = false,
#                 framestyle = :box,
#                 yticks = ([0.05, 0.1, 0.15, 0.2], ["0.05", "0.10", "0.15", "0.20"]),
#                 xticks = (-10:5:40, ["-10", "", "0", "", "10", "", "20", "", "30", "", "40"])
#                 )
#     plt6 = plot(x, evo_b[inds[3]+1].data[y_inds],
#                 title = "\$b_{$(inds[3]+1)}\$",
#                 xlims = (-10, 45),
#                 ylims = (0,0.25),
#                 legend = false,
#                 framestyle = :box,
#                 yticks = ([0.05, 0.1, 0.15, 0.2], ["0.05", "0.10", "0.15", "0.20"]),
#                 xticks = (-10:5:40, ["-10", "", "0", "", "10", "", "20", "", "30", "", "40"])
#                 )
#     plt7 = plot(x, evo_a[inds[4]].data[y_inds],
#                 title = "\$a_{$(inds[4])}\$",
#                 xlims = (-10, 45),
#                 ylims = (0,0.25),
#                 legend = false,
#                 framestyle = :box,
#                 yticks = ([0.05, 0.1, 0.15, 0.2], ["0.05", "0.10", "0.15", "0.20"]),
#                 xticks = (-10:5:40, ["-10", "", "0", "", "10", "", "20", "", "30", "", "40"])
#                 )
#     plt8 = plot(x, evo_b[inds[4]+1].data[y_inds],
#                 title = "\$b_{$(inds[4]+1)}\$",
#                 xlims = (-10, 45),
#                 ylims = (0,0.25),
#                 legend = false,
#                 framestyle = :box,
#                 yticks = ([0.05, 0.1, 0.15, 0.2], ["0.05", "0.10", "0.15", "0.20"]),
#                 xticks = (-10:5:40, ["-10", "", "0", "", "10", "", "20", "", "30", "", "40"])
#                 )
#     plt9 = plot(x, evo_a[inds[5]].data[y_inds],
#                 title = "\$a_{$(inds[5])}\$",
#                 xlims = (-10, 45),
#                 ylims = (0,0.25),
#                 legend = false,
#                 framestyle = :box,
#                 yticks = ([0.05, 0.1, 0.15, 0.2], ["0.05", "0.10", "0.15", "0.20"]),
#                 xticks = (-10:5:40, ["-10", "", "0", "", "10", "", "20", "", "30", "", "40"])
#                 )
#     plt10 = plot(x, evo_b[inds[5]+1].data[y_inds],
#                 title = "\$b_{$(inds[5]+1)}\$",
#                 xlims = (-10, 45),
#                 ylims = (0,0.25),
#                 legend = false,
#                 framestyle = :box,
#                 yticks = ([0.05, 0.1, 0.15, 0.2], ["0.05", "0.10", "0.15", "0.20"]),
#                 xticks = (-10:5:40, ["-10", "", "0", "", "10", "", "20", "", "30", "", "40"])
#                 )
#     plt = plot(plt1, plt2, plt3, plt4, plt5, plt6, plt7, plt8, plt9, plt10, layout = (5,2), size = (600, 900))

#     return plt, evo_a, evo_b
# end

# #########################################################################
# ########## below is old code, delete after verifying the above ##########
# #########################################################################

# # # @assert l_max > 1
# # # @assert 0 < ε < 1
# # # channel ∈ (:BEC, :BSC, :BIAWGN) || throw(ArgumentError("Channel not yet implemented"))
# # # decoder ∈ (:BEC, :A, :SP) || throw(ArgumentError("Decoder not supported"))
# # # if channel == :BEC
# # #     decoder == :BEC || throw(ArgumentError("The only decoder supported for the BEC channel is :BEC"))
# # # end

# # # R, x = PolynomialRing(RealField(), :x)

# # function optimal_lambda_and_rho(l_max::Int, r_max::Int, real_param::Float64, var_type::Symbol)
# #     var_type ∈ (:r, :ε) || throw(ArgumentError("var_type must be :r for target rate or :ε for threshold"))
# #     var_type == :r && real_param >= 1 - 2/r_max && throw(ArgumentError("This rate is unachieveable with the given r_max."))
# #     # TODO: check for when var_type == :ε as well

# #     tolerance = 1e-9

# #     # initial guess: ρ(x) = x^(r_max - 1)
# #     ρ = zeros(r_max); ρ[end] = 1
# #     # this makes sense to me, but we could choose for some c ∈ (0, 1]
# #     #   ρ(x) = (1 - c) * x^(r_max - 2) + c * x^(r_max - 1)
# #     # which would be given by the code:
# #     # ρ = zeros(r_max); c = 0.5; ρ[end - 1] = 1 - c; ρ[end] = c;

# #     # if we need an initial λ, this would be it:
# #     # λ, _, _ = _optimal_distributions(ρ, :ρ, l_max, real_param, var_type)

# #     # solve until convergence for each ε, see if rates match, else change ε and repeat
# #     if var_type == :r
# #         max_iters = 100
# #         high = 1.0
# #         low = 0.0
# #         mid = 0.5

# #         countinner = 0
# #         λ, _ = _find_lambda_given_rho(ρ, mid, l_max)
# #         ρ, _ = _find_rho_given_lambda(λ, mid, r_max)
# #         λprev = copy(λ)
# #         ρprev = copy(ρ)
# #         convergedinner = false
# #         while countinner <= max_iters
# #             countinner += 1
# #             λ, _ = _find_lambda_given_rho(ρ, mid, l_max)
# #             ρ, _ = _find_rho_given_lambda(λ, mid, r_max)
# #             normλ = _L2_dist_sq(λ, λprev)
# #             normρ = _L2_dist_sq(ρ, ρprev)
# #             normλ <= tolerance && normρ <= tolerance && (convergedinner = true; break;)
# #             λprev .= λ
# #             ρprev .= ρ
# #         end
# #         if !convergedinner
# #             # TODO: better error here, just putting something for now
# #             error("inner convergence failed")
# #         end
# #         sol_rate = 1 - _integrate_poly_0_1(ρ) / _integrate_poly_0_1(λ)
# #         Δ = sol_rate - real_param

# #         converged = abs(Δ) <= tolerance
# #         count = 0
# #         while count <= max_iters && !converged
# #             count += 1
# #             Δ > 0 ? (low = mid;) : (high = mid;)
# #             mid = (high + low) / 2

# #             countinner = 0
# #             convergedinner = false
# #             while countinner <= max_iters
# #                 countinner += 1
# #                 λ, _ = _find_lambda_given_rho(ρ, mid, l_max)
# #                 ρ, _ = _find_rho_given_lambda(λ, mid, r_max)
# #                 normλ = _L2_dist_sq(λ, λprev)
# #                 normρ = _L2_dist_sq(ρ, ρprev)
# #                 normλ <= tolerance && normρ <= tolerance && (convergedinner = true; break;)
# #                 λprev .= λ
# #                 ρprev .= ρ
# #             end
# #             if !convergedinner
# #                 # TODO: better error here, just putting something for now
# #                 error("inner convergence failed")
# #             end
# #             converged = abs(Δ) <= tolerance
# #         end
# #         # TODO: better error here, just putting something for now
# #         converged ? (return λ, ρ, mid;) : error("outer convergence failed")
# #     else # var_type == :ε
# #         countinner = 0
# #         λ, _ = _find_lambda_given_rho(ρ, real_param, l_max)
# #         ρ, _ = _find_rho_given_lambda(λ, real_param, r_max)
# #         λprev = copy(λ)
# #         ρprev = copy(ρ)
# #         convergedinner = false
# #         while countinner <= max_iters
# #             countinner += 1
# #             λ, _ = _find_lambda_given_rho(ρ, real_param, l_max)
# #             ρ, _ = _find_rho_given_lambda(λ, real_param, r_max)
# #             normλ = _L2_dist_sq(λ, λprev)
# #             normρ = _L2_dist_sq(ρ, ρprev)
# #             normλ <= tolerance && normρ <= tolerance && (convergedinner = true; break;)
# #             λprev .= λ
# #             ρprev .= ρ
# #         end
# #         convergedinner && (return λ, ρ, 1 - _integrate_poly_0_1(ρ) / _integrate_poly_0_1(λ))
# #     end
# # end

# # function optimal_threshold(λ, ρ)
# #     high = 1.0
# #     low = 0.0
# #     Δ = 1
# #     tolerance = 1e-9
# #     while Δ > tolerance
# #         mid = (high + low) / 2
# #         # some way to evaluate f(x) on this range
# #         # for x in 0.001:001:1
# #         #     # l_max
# #         #     mid * sum(λ[i] * (1 - _poly_eval(1 - x, ρ_vec))^i for i in 1:l_max - 1) - x
# #         # end

# #         # some loop
# #         f(low) * f(high) < 0 ? (high = mid;) : (low = mid;)
# #         Δ = abs(f(mid))
# #     end
# #     return mid
# # end
