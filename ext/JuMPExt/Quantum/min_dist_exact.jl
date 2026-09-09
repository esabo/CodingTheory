# Copyright (c) 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
    _minimum_distance_css_ILP(H, logical_checks; max_d, verbose, time_limit_sec)

Exactly minimize the weight of a binary vector in `ker(H)` whose logical
label under `logical_checks` is nonzero.
"""
function _minimum_distance_css_ILP(
    H, logical_checks;
    max_d::Int=size(H, 2), verbose::Bool=false,
    time_limit_sec::Float64=300.0
)
    size(H, 2) == size(logical_checks, 2) ||
        throw(ArgumentError("The parity-check and logical-check matrices must have the same number of columns."))
    0 <= max_d <= size(H, 2) ||
        throw(DomainError(max_d, "`max_d` must lie between zero and the code length."))
    time_limit_sec > 0 ||
        throw(DomainError(time_limit_sec, "`time_limit_sec` must be positive."))

    H_int = CodingTheory._convert_binary_to_int_matrix(H)
    L_int = CodingTheory._convert_binary_to_int_matrix(logical_checks)
    m, n = size(H_int)
    ell = size(L_int, 1)
    iszero(ell) && return -1, zeros(Int, n)

    model = Model(GLPK.Optimizer)
    verbose || set_silent(model)
    set_time_limit_sec(model, time_limit_sec)

    @variable(model, x[1:n], Bin)
    @variable(model, syndrome_quotient[1:m] >= 0, Int)
    @variable(model, logical_label[1:ell], Bin)
    @variable(model, logical_quotient[1:ell] >= 0, Int)

    for r in 1:m
        @constraint(model,
            sum(H_int[r, c] * x[c] for c in 1:n) ==
            2 * syndrome_quotient[r])
    end
    for r in 1:ell
        @constraint(model,
            sum(L_int[r, c] * x[c] for c in 1:n) ==
            2 * logical_quotient[r] + logical_label[r])
    end
    @constraint(model, sum(logical_label) >= 1)
    @constraint(model, sum(x) <= max_d)
    @objective(model, Min, sum(x))

    optimize!(model)
    status = termination_status(model)
    if status == JuMP.MOI.INFEASIBLE
        return -1, zeros(Int, n)
    elseif status != JuMP.MOI.OPTIMAL
        error("CSS minimum-distance ILP did not finish exactly (termination status: $status).")
    end

    witness = round.(Int, value.(x))
    return sum(witness), witness
end
