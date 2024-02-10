# Copyright (c) 2023 - 2024 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

function _find_lambda_given_rho(ρ::Union{Vector{Float64}, CodingTheory.Oscar.PolyRingElem},
    ε::Float64, l_max::Int; Δ = 0.001)

    model = Model(GLPK.Optimizer)
    @variable(model, λ[1:l_max - 1] >= 0)
    @constraint(model, sum(λ) == 1)
    for x in 0:Δ:1
        @constraint(model, ε * sum(λ[i] * (1 - CodingTheory._poly_eval(1 - x, ρ))^i for i in
            1:l_max - 1) - x <= 0)
    end
    @constraint(model, ε * CodingTheory._d_poly_eval(1, ρ) * λ[1] <= 1)
    @objective(model, Max, sum(λ[i] / (i + 1) for i in 1:l_max - 1))
    optimize!(model)
    termination_status(model) == MOI.INFEASIBLE && throw(DomainError("No solution exists"))
    @assert termination_status(model) == MOI.OPTIMAL "Didn't find an optimal point"
    return [0.0; value.(λ)], model
end

function _find_rho_given_lambda(λ::Union{Vector{Float64}, CodingTheory.Oscar.PolyRingElem},
    ε::Float64, r_max::Int; Δ = 0.001)

    model = Model(GLPK.Optimizer)
    @variable(model, ρ[1:r_max - 1] >= 0)
    @constraint(model, sum(ρ) == 1)
    for x in 0:Δ:1
        @constraint(model, 1 - x - sum(ρ[i] * (1 - ε * CodingTheory._poly_eval(x, λ))^i for i in
            1:r_max - 1) <= 0)
    end
    temp = isa(λ, PolyRingElem) ? Float64(coeff(λ, 1)) : λ[2]
    @constraint(model, ε * sum(i * ρ[i] for i in eachindex(ρ)) * temp <= 1)
    @objective(model, Min, sum(ρ[i] / (i + 1) for i in 1:r_max - 1))
    optimize!(model)
    termination_status(model) == MOI.INFEASIBLE && throw(DomainError("No solution exists"))
    @assert termination_status(model) == MOI.OPTIMAL "Didn't find an optimal point"
    return [0.0; value.(ρ)], model
end

"""
    optimal_lambda(ρ, l_max, param, var_type; Δ = 1e-3)

Find the optimal variable node distribution `λ` given the check node
distribution `ρ`, maximum variable node degree `l_max`, and target parameter
`param` which refers to threshold if `var_type == :ε` or rate if `var_type == :r`.

# Notes
* `Δ` refers to the step size for `x` in the LP to solve for `λ`.
"""
function CodingTheory.optimal_lambda(ρ, l_max::Int, param::Float64, var_type::Symbol; Δ = 0.001)
    var_type ∈ (:r, :ε) || throw(ArgumentError("var_type must be :r for target rate or :ε for threshold"))
    var_type == :r && param >= 1 - 2 * CodingTheory._integrate_poly_0_1(ρ) && throw(
        ArgumentError("This rate is unachieveable with the given ρ."))
    # TODO: check for when var_type == :ε as well

    λ_vec, r, ε =  _optimal_distributions(ρ, :ρ, l_max, param, var_type, Δλ = 0.001)
    _, x = CodingTheory.Oscar.PolynomialRing(CodingTheory.Oscar.RealField(), :x)
    λ = sum(c * x^(i - 1) for (i, c) in enumerate(λ_vec))
    return (λ = λ, r = r, ε = ε)
end

"""
    optimal_rho(λ, r_max, param, var_type; Δ = 1e-3)

Find the optimal check node distribution `ρ` given the variable node
distribution `λ`, maximum check node degree `r_max`, and target parameter `param`
which refers to threshold if `var_type == :ε` or rate if `var_type == :r`.

# Notes
* `Δ` refers to the step size for x in the LP to solve for ρ.
"""
function CodingTheory.optimal_rho(λ, r_max::Int, param::Float64, var_type::Symbol, Δ = 1e-3)
    var_type ∈ (:r, :ε) || throw(ArgumentError("var_type must be :r for target rate or :ε for threshold"))
    var_type == :r && param >= 1 - 1 / (r_max * CodingTheory._integrate_poly_0_1(λ)) &&
        throw(ArgumentError("This rate is unachieveable with the given r_max and λ."))
    # TODO: check for when var_type == :ε as well

    ρ_vec, r, ε = _optimal_distributions(λ, :λ, r_max, param, var_type, Δρ = 0.001)
    _, x = CodingTheory.Oscar.PolynomialRing(CodingTheory.Oscar.RealField(), :x)
    ρ = sum(c * x^(i - 1) for (i, c) in enumerate(ρ_vec))
    return (ρ = ρ, r = r, ε = ε)
end

# TODO: add a type on poly
function _optimal_distributions(poly, poly_type::Symbol, var_max::Int, real_param::Float64,
    var_type::Symbol; Δλ = 0.001, Δρ = 0.001)

    int_poly = CodingTheory._integrate_poly_0_1(poly)
    if var_type == :r
        tolerance = 1e-6
        max_iters = 100
        high = 1.0
        low = 0.0
        mid = 0.5
        count = 0
        Δ = -Inf
        while !(0 <= Δ <= tolerance) && count < max_iters
            count += 1
            mid = (high + low) / 2
            Δ, sol = try
                if poly_type == :ρ
                    sol, _ = _find_lambda_given_rho(poly, mid, var_max, Δ = Δλ)
                    sol_rate = 1 - int_poly / CodingTheory._integrate_poly_0_1(sol)
                else
                    sol, _ = _find_rho_given_lambda(poly, mid, var_max, Δ = Δρ)
                    sol_rate = 1 - CodingTheory._integrate_poly_0_1(sol) / int_poly
                end
                (sol_rate - real_param, sol)
            catch
                (-Inf, Float64[])
            end
            Δ > 0 ? (low = mid;) : (high = mid;)
        end
        0 <= Δ <= tolerance || error("Solution for $(poly_type == :ρ ? :λ : :ρ) did not converge in $max_iters iterations")
        return sol, sol_rate, mid
    else
        if poly_type == :ρ
            sol, _ = _find_lambda_given_rho(poly, real_param, var_max, Δ = Δλ)
            sol_rate = 1 - int_poly / CodingTheory._integrate_poly_0_1(sol)
        else
            sol, _ = _find_rho_given_lambda(poly, real_param, var_max, Δ = Δρ)
            sol_rate = 1 - CodingTheory._integrate_poly_0_1(sol) / int_poly
        end
        return sol, sol_rate, real_param
    end
end

function _find_lambda_and_rho(l_max::Int, r_max::Int, ε::Float64; Δρ = 0.01, Δλ = 0.001)
    # TODO: probably try a smallish Δρ (0.05?) and then bisect in the interval where it's best.
    c = 0:Δρ:1
    iters = length(c)
    ρ = zeros(r_max, iters)
    λ = zeros(l_max, iters)
    rates = fill(-Inf, iters)
    for i in 1:iters
        ρ[end, i] = c[i]
        ρ[end - 1, i] = (1 - c[i])
        try
            λ[:, i] .= _find_lambda_given_rho(ρ[:, i], ε, l_max, Δ = Δλ)[1]
            rates[i] = 1 - CodingTheory._integrate_poly_0_1(ρ[:, i]) /
                CodingTheory._integrate_poly_0_1(λ[:, i])
        catch end
    end
    i = argmax(rates)
    isfinite(rates[i]) || throw(ArgumentError("No solution for given parameters"))
    return λ[:, i], ρ[:, i], rates[i]
end

"""
    optimal_lambda_and_rho(l_max, r_max, param, var_type; Δρ = 1e-2, Δλ = 1e-3)

Find the optimal distribution pair λ, ρ given the `param`, where `param` is
either a threshold if `var_type == :ε` or a target rate if `var_type == :r`.

# Notes
* `Δρ` gives the step size for possible values of `c` where
    ρ = (1 - c) * x^(r_max - 2) + c * x^(r_max - 1)
* `Δλ` gives the step size for values of `x` in the LP for finding λ given ρ.
"""
function CodingTheory.optimal_lambda_and_rho(l_max::Int, r_max::Int, param::Float64,
    var_type::Symbol; Δρ = 0.01, Δλ = 0.001)

    if var_type == :r
        tolerance = 1e-6
        max_iters = 100
        high = 1.0
        low = 0.0
        mid = 0.5
        count = 0
        Δ = -Inf
        while !(0 <= Δ <= tolerance) && count < max_iters
            count += 1
            mid = (high + low) / 2
            Δ, λ_vec, ρ_vec = try
                λ, ρ, sol_rate = _find_lambda_and_rho(l_max, r_max, mid, Δρ = Δρ, Δλ = Δλ)
                (sol_rate - param, λ, ρ)
            catch
                (-Inf, Float64[], Float64[])
            end
            Δ > 0 ? (low = mid;) : (high = mid;)
        end
        0 <= Δ <= tolerance || error("Solution for $(poly_type == :ρ ? :λ : :ρ) did not converge in $max_iters iterations")
        _, x = CodingTheory.Oscar.PolynomialRing(CodingTheory.Oscar.RealField(), :x)
        λ = sum(c * x^(i - 1) for (i, c) in enumerate(λ_vec))
        ρ = sum(c * x^(i - 1) for (i, c) in enumerate(ρ_vec))
        return (λ = λ, ρ = ρ, r = sol_rate, ε = mid)
    elseif var_type == :ε
        λ_vec, ρ_vec, sol_rate = _find_lambda_and_rho(l_max, r_max, param, Δρ = Δρ, Δλ = Δλ)
        _, x = CodingTheory.Oscar.PolynomialRing(CodingTheory.Oscar.RealField(), :x)
        λ = sum(c * x^(i - 1) for (i, c) in enumerate(λ_vec))
        ρ = sum(c * x^(i - 1) for (i, c) in enumerate(ρ_vec))
        return (λ = λ, ρ = ρ, r = sol_rate, ε = param)
    else
        throw(ArgumentError("var_type must be :r or :ε"))
    end
end
