module QuantumBoundsExt

using CodingTheory
using JuMP
using Tulip

import CodingTheory:
    QuantumLPResult,
    _quantum_CSS_dimension_LP_bound,
    _quantum_CSS_weight_enumerator_LP,
    _quantum_stabilizer_dimension_LP_bound,
    _quantum_stabilizer_generator_weight_LP_bound,
    _quantum_weight_enumerator_LP,
    quantum_Krawtchouk_matrix,
    quantum_stabilizer_generator_weight_lower_bound

const MOI = JuMP.MOI

struct ExactRow
    coefficients::Vector{BigInt}
    sense::Symbol
    rhs::BigInt
end

function _row(
    n::Int, terms, sense::Symbol, rhs::Integer=0,
)
    coefficients = zeros(BigInt, n + 1)
    for (index, value) in terms
        coefficients[index + 1] += value
    end
    return ExactRow(coefficients, sense, BigInt(rhs))
end

function _exact_rows(
    n::Int, k::Int, d::Int;
    formulation::Symbol=:standard,
    check_weight::Union{Nothing, Int}=nothing,
    num_max_weight_generators::Union{Nothing, Int}=nothing,
    parity::Symbol=:auto,
    connected::Bool=false,
    cumulative_lower_bounds=nothing,
    include_shadow::Bool=true,
)
    formulation in (:standard, :coarse, :refined) ||
        throw(ArgumentError(
            "Expected formulation=:standard, :coarse, or :refined."))
    r = n - k
    rows = ExactRow[]
    push!(rows, _row(n, [(0, 1)], :eq, 1))

    M = quantum_Krawtchouk_matrix(n)
    signed_M = include_shadow ?
        quantum_Krawtchouk_matrix(n; signed_columns=true) : nothing
    stabilizer_size = BigInt(2)^r
    for j in 0:n
        coefficients = copy(M[j + 1, :])
        coefficients[j + 1] -= stabilizer_size
        push!(rows, ExactRow(
            coefficients, j < d ? :eq : :ge, BigInt(0)))
        include_shadow &&
            push!(rows, ExactRow(copy(signed_M[j + 1, :]), :ge, BigInt(0)))
    end
    if !isnothing(cumulative_lower_bounds)
        for (cutoff, count) in cumulative_lower_bounds
            0 <= cutoff <= n ||
                throw(DomainError(cutoff, "A cumulative cutoff must lie in 0:n."))
            count >= 0 ||
                throw(DomainError(count, "A cumulative lower bound cannot be negative."))
            push!(rows, _row(
                n, ((i, 1) for i in 0:Int(cutoff)), :ge, count))
        end
    end

    formulation == :standard && return rows
    isnothing(check_weight) &&
        throw(ArgumentError("check_weight is required for $formulation."))
    w = check_weight
    1 <= w <= n ||
        throw(DomainError(w, "check_weight must lie in 1:n."))

    for C in 1:fld(n, w)
        rhs = sum(binomial(BigInt(r), c) for c in 1:C; init=BigInt(0))
        push!(rows, _row(
            n, ((i, 1) for i in 1:min(C * w, n)), :ge, rhs))
    end
    formulation == :coarse && return rows

    isnothing(num_max_weight_generators) &&
        throw(ArgumentError(
            "num_max_weight_generators is required for formulation=:refined."))
    y = num_max_weight_generators
    1 <= y <= r ||
        throw(DomainError(y, "num_max_weight_generators must lie in 1:(n-k)."))
    push!(rows, _row(n, [(w, 1)], :ge, y))
    push!(rows, _row(
        n, ((i, 1) for i in 0:(w - 1)), :le, BigInt(2)^(r - y)))

    weight_buckets = zeros(BigInt, n + 1)
    for p in 0:(r - y), q in 0:y
        weight = p * (w - 1) + q * w
        weight <= n || continue
        weight_buckets[weight + 1] +=
            binomial(BigInt(r - y), p) * binomial(BigInt(y), q)
    end
    cumulative_counts = cumsum(weight_buckets)
    for cutoff in (w - 1):(n - 1)
        rhs = cumulative_counts[cutoff + 1]
        push!(rows, _row(
            n, ((i, 1) for i in 0:cutoff), :ge, rhs))
    end

    parity in (:auto, :half, :all_even) ||
        throw(ArgumentError("Expected parity=:auto, :half, or :all_even."))
    parity == :auto && iseven(w) && y < r &&
        throw(ArgumentError(
            "Even check weight with y < n-k requires testing both " *
            "parity=:half and parity=:all_even."))
    resolved_parity = parity == :auto ?
        (isodd(w) ? :half : :all_even) : parity
    rhs = resolved_parity == :half ? BigInt(2)^(r - 1) : BigInt(2)^r
    push!(rows, _row(
        n, ((i, 1) for i in 0:2:n), :eq, rhs))

    if k < n - 1
        shorter_bound =
            quantum_stabilizer_generator_weight_lower_bound(n - 1, k)
        w < shorter_bound &&
            push!(rows, _row(n, [(1, 1)], :eq, 0))
    end
    connected && push!(rows, _row(
        n, ((i, 1) for i in 1:min(2w - 2, n)), :ge, 2r - 1))
    return rows
end

function _scaled(row::ExactRow)
    scale = max(abs(row.rhs), maximum(abs, row.coefficients))
    iszero(scale) && return row.coefficients, row.rhs
    return row.coefficients .// scale, row.rhs // scale
end

function _normalized_violation(
    rows::Vector{ExactRow}, solution::Vector{BigFloat},
)
    violation = max(BigFloat(0), maximum(max(-x, 0) for x in solution))
    for row in rows
        coefficients, rhs = _scaled(row)
        residual = sum(BigFloat(coefficients[i]) * solution[i]
            for i in eachindex(solution)) - BigFloat(rhs)
        row_violation = row.sense == :eq ? abs(residual) :
            row.sense == :ge ? max(-residual, 0) : max(residual, 0)
        violation = max(violation, row_violation)
    end
    return violation
end

function _solve_once(
    rows::Vector{ExactRow}, num_variables::Int, precision::Int;
    model_hook=nothing,
)
    return setprecision(BigFloat, precision) do
        model = GenericModel{BigFloat}(Tulip.Optimizer{BigFloat})
        set_silent(model)
        @variable(model, A[1:num_variables] >= 0)
        for row in rows
            coefficients, rhs = _scaled(row)
            expression = sum(
                BigFloat(coefficients[i]) * A[i] for i in 1:num_variables)
            rhs_float = BigFloat(rhs)
            if row.sense == :eq
                @constraint(model, expression == rhs_float)
            elseif row.sense == :ge
                @constraint(model, expression >= rhs_float)
            else
                @constraint(model, expression <= rhs_float)
            end
        end
        isnothing(model_hook) || model_hook(model, A)
        @objective(model, Min, zero(BigFloat))
        optimize!(model)
        status = termination_status(model)
        if status in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL)
            solution = BigFloat[value(A[i]) for i in 1:num_variables]
            violation = _normalized_violation(rows, solution)
            tolerance = BigFloat(2)^(-fld(precision, 3))
            verified_status = violation <= tolerance ? :feasible : :unknown
            return verified_status, Symbol(string(status)), solution, violation
        elseif status in (MOI.INFEASIBLE, MOI.ALMOST_INFEASIBLE)
            return :infeasible_numerical, Symbol(string(status)),
                nothing, BigFloat(Inf)
        end
        return :unknown, Symbol(string(status)), nothing, BigFloat(Inf)
    end
end

function _CSS_exact_rows(
    n::Int, k_X::Int, k_Z::Int, d::Int;
    check_weight::Union{Nothing, Int}=nothing,
    exclude_weight_one::Bool=false,
)
    num_variables = 2(n + 1)
    x(i) = i + 1
    z(i) = n + 2 + i
    row(terms, sense, rhs=0) = begin
        coefficients = zeros(BigInt, num_variables)
        for (index, value) in terms
            coefficients[index] += value
        end
        ExactRow(coefficients, sense, BigInt(rhs))
    end

    rows = ExactRow[
        row([(x(0), 1)], :eq, 1),
        row([(z(0), 1)], :eq, 1),
        row(((x(i), 1) for i in 0:n), :eq, BigInt(2)^(n - k_X)),
        row(((z(i), 1) for i in 0:n), :eq, BigInt(2)^(n - k_Z)),
    ]
    M = quantum_Krawtchouk_matrix(n; alphabet_size=2)
    size_X = BigInt(2)^(n - k_X)
    size_Z = BigInt(2)^(n - k_Z)
    for ell in 0:n
        # B^Z_ell - A^X_ell and B^X_ell - A^Z_ell, with the
        # MacWilliams denominators cleared exactly.
        z_to_x = zeros(BigInt, num_variables)
        x_to_z = zeros(BigInt, num_variables)
        for j in 0:n
            z_to_x[z(j)] = M[ell + 1, j + 1]
            x_to_z[x(j)] = M[ell + 1, j + 1]
        end
        z_to_x[x(ell)] -= size_Z
        x_to_z[z(ell)] -= size_X
        sense = 1 <= ell < d ? :eq : :ge
        push!(rows, ExactRow(z_to_x, sense, BigInt(0)))
        push!(rows, ExactRow(x_to_z, sense, BigInt(0)))
    end

    if !isnothing(check_weight)
        w = check_weight
        1 <= w <= n ||
            throw(DomainError(w, "check_weight must lie in 1:n."))
        for m in 0:fld(n, w)
            rhs_X = sum(
                binomial(BigInt(n - k_X), j) for j in 0:m;
                init=BigInt(0))
            rhs_Z = sum(
                binomial(BigInt(n - k_Z), j) for j in 0:m;
                init=BigInt(0))
            push!(rows, row(
                ((x(i), 1) for i in 0:min(m * w, n)), :ge, rhs_X))
            push!(rows, row(
                ((z(i), 1) for i in 0:min(m * w, n)), :ge, rhs_Z))
        end
    end
    if exclude_weight_one
        push!(rows, row([(x(1), 1)], :eq, 0))
        push!(rows, row([(z(1), 1)], :eq, 0))
    end
    return rows
end

function _quantum_CSS_weight_enumerator_LP(
    n::Int, k_X::Int, k_Z::Int, d::Int;
    check_weight::Union{Nothing, Integer}=nothing,
    exclude_weight_one::Bool=false,
    precisions=(256, 512),
    model_hook=nothing,
)
    rows = _CSS_exact_rows(
        n, k_X, k_Z, d;
        check_weight=isnothing(check_weight) ? nothing : Int(check_weight),
        exclude_weight_one=exclude_weight_one,
    )
    result = nothing
    for precision in precisions
        precision >= 64 ||
            throw(DomainError(precision, "LP precision must be at least 64 bits."))
        status, termination, enumerator, violation =
            _solve_once(
                rows, 2(n + 1), Int(precision); model_hook=model_hook)
        result = QuantumLPResult(
            status, :css, Int(precision), termination, enumerator, violation)
        status == :feasible && return result
    end
    return result
end

function _quantum_weight_enumerator_LP(
    n::Int, k::Int, d::Int;
    formulation::Symbol=:standard,
    check_weight::Union{Nothing, Integer}=nothing,
    num_max_weight_generators::Union{Nothing, Integer}=nothing,
    parity::Symbol=:auto,
    connected::Bool=false,
    cumulative_lower_bounds=nothing,
    include_shadow::Bool=true,
    precisions=(256, 512),
    model_hook=nothing,
)
    rows = _exact_rows(
        n, k, d;
        formulation=formulation,
        check_weight=isnothing(check_weight) ? nothing : Int(check_weight),
        num_max_weight_generators=isnothing(num_max_weight_generators) ?
            nothing : Int(num_max_weight_generators),
        parity=parity,
        connected=connected,
        cumulative_lower_bounds=cumulative_lower_bounds,
        include_shadow=include_shadow,
    )
    attempts = QuantumLPResult[]
    for precision in precisions
        precision >= 64 ||
            throw(DomainError(precision, "LP precision must be at least 64 bits."))
        status, termination, enumerator, violation =
            _solve_once(rows, n + 1, Int(precision); model_hook=model_hook)
        push!(attempts, QuantumLPResult(
            status, formulation, Int(precision), termination,
            enumerator, violation))
        status == :feasible && return attempts[end]
        status == :unknown && continue
    end
    return attempts[end]
end

function _quantum_stabilizer_generator_weight_LP_bound(
    n::Int, k::Int, d::Int;
    max_weight::Integer=n,
    connected::Bool=false,
    precisions=(256, 512),
)
    r = n - k
    first_weight =
        Int(quantum_stabilizer_generator_weight_lower_bound(n, k))
    trials = NamedTuple[]
    for w in first_weight:min(Int(max_weight), n)
        minimum_y = max(1, 2n - (w - 1) * r)
        unknown = false
        for y in minimum_y:r
            parities = iseven(w) && y < r ? (:half, :all_even) : (:auto,)
            for parity in parities
                result = _quantum_weight_enumerator_LP(
                    n, k, d;
                    formulation=:refined,
                    check_weight=w,
                    num_max_weight_generators=y,
                    parity=parity,
                    connected=connected,
                    precisions=precisions,
                )
                push!(trials, (; check_weight=w, y, parity, result))
                result.status == :feasible &&
                    return (; lower_bound=w, status=:not_excluded, trials)
                unknown |= result.status == :unknown
            end
        end
        unknown &&
            return (; lower_bound=w, status=:unknown, trials)
    end
    return (;
        lower_bound=min(Int(max_weight), n) + 1,
        status=:infeasible_numerical,
        trials,
    )
end

function _quantum_stabilizer_dimension_LP_bound(
    n::Int, d::Int, check_weight::Int;
    exclude_weight_one::Bool=false,
    include_shadow::Bool=false,
    precisions=(256, 512),
)
    max_k = min(n - 1, n - 2d + 2)
    trials = NamedTuple[]
    for k in max_k:-1:1
        hook = exclude_weight_one ?
            ((model, A) -> @constraint(model, A[2] == 0)) : nothing
        result = _quantum_weight_enumerator_LP(
            n, k, d;
            formulation=:coarse,
            check_weight=check_weight,
            include_shadow=include_shadow,
            precisions=precisions,
            model_hook=hook,
        )
        push!(trials, (; k, result))
        result.status == :feasible &&
            return (; upper_bound=k, status=:not_excluded_boundary, trials)
        result.status == :unknown &&
            return (; upper_bound=k, status=:unknown, trials)
    end
    return (; upper_bound=0, status=:infeasible_numerical, trials)
end

function _quantum_CSS_dimension_LP_bound(
    n::Int, d::Int, check_weight::Int;
    exclude_weight_one::Bool=false,
    precisions=(256, 512),
)
    max_k = min(n - 1, n - 2d + 2)
    trials = NamedTuple[]
    for k in max_k:-1:1
        unknown = false
        for k_Z in max(0, k):fld(n + k, 2)
            k_X = n + k - k_Z
            result = _quantum_CSS_weight_enumerator_LP(
                n, k_X, k_Z, d;
                check_weight=check_weight,
                exclude_weight_one=exclude_weight_one,
                precisions=precisions,
            )
            push!(trials, (; k, k_X, k_Z, result))
            result.status == :feasible &&
                return (; upper_bound=k, status=:not_excluded_boundary, trials)
            unknown |= result.status == :unknown
        end
        unknown &&
            return (; upper_bound=k, status=:unknown, trials)
    end
    return (; upper_bound=0, status=:infeasible_numerical, trials)
end

end
