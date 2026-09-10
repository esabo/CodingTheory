function _solve_css_minimum_distance_ilp(
    optimizer_factory, H, logical_checks;
    min_d::Int=1, max_d::Int=size(H, 2), verbose::Bool=false,
    time_limit_sec::Union{Nothing, Float64}=nothing,
    parity_cut_max_degree::Int=10, verbose_attributes=Pair[],
    threads::Int=0, cyclic_period::Union{Nothing, Int}=nothing
)
    size(H, 2) == size(logical_checks, 2) ||
        throw(ArgumentError("The parity-check and logical-check matrices must have the same number of columns."))
    0 <= max_d <= size(H, 2) ||
        throw(DomainError(max_d, "`max_d` must lie between zero and the code length."))
    1 <= min_d <= max_d + 1 ||
        throw(DomainError(min_d, "`min_d` must lie between one and `max_d + 1`."))
    isnothing(time_limit_sec) || time_limit_sec > 0 ||
        throw(DomainError(time_limit_sec, "`time_limit_sec` must be positive."))
    0 <= parity_cut_max_degree <= 20 ||
        throw(DomainError(parity_cut_max_degree,
            "`parity_cut_max_degree` must lie between zero and 20."))
    threads >= 0 ||
        throw(DomainError(threads, "`threads` must be nonnegative."))

    H_int = CodingTheory._convert_binary_to_int_matrix(H)
    L_int = CodingTheory._convert_binary_to_int_matrix(logical_checks)
    m, n = size(H_int)
    ell = size(L_int, 1)
    iszero(ell) && return -1, zeros(Int, n), :infeasible

    model = Model(optimizer_factory)
    if verbose
        for (attribute, value) in verbose_attributes
            set_attribute(model, attribute, value)
        end
    else
        set_silent(model)
    end
    isnothing(time_limit_sec) || set_time_limit_sec(model, time_limit_sec)
    if threads > 0
        set_attribute(model, "threads", threads)
        set_attribute(model, "parallel", "on")
    end
    try
        set_attribute(model, "mip_detect_symmetry", true)
    catch
        # Older HiGHS builds do not expose this option.
    end

    @variable(model, x[1:n], Bin)
    @variable(model, syndrome_quotient[1:m] >= 0, Int)
    @variable(model, logical_label[1:ell], Bin)
    @variable(model, logical_quotient[1:ell] >= 0, Int)

    parity_cut_count = 0
    for r in 1:m
        support = findall(!iszero, @view H_int[r, :])
        @constraint(model,
            sum(x[c] for c in support) == 2 * syndrome_quotient[r])
        @constraint(model, syndrome_quotient[r] <= length(support) ÷ 2)

        degree = length(support)
        if 1 <= degree <= parity_cut_max_degree
            for mask in UInt(1):((UInt(1) << degree) - UInt(1))
                isodd(count_ones(mask)) || continue
                selected = count_ones(mask)
                @constraint(model,
                    sum(((mask >> (j - 1)) & UInt(1) == UInt(1) ? 1 : -1) *
                        x[support[j]] for j in 1:degree) <= selected - 1)
                parity_cut_count += 1
            end
        end
    end
    for r in 1:ell
        support = findall(!iszero, @view L_int[r, :])
        @constraint(model,
            sum(x[c] for c in support) ==
            2 * logical_quotient[r] + logical_label[r])
        @constraint(model, logical_quotient[r] <= length(support) ÷ 2)
    end
    @constraint(model, sum(logical_label) >= 1)
    @constraint(model, sum(x) >= min_d)
    @constraint(model, sum(x) <= max_d)
    if cyclic_period !== nothing
        period = Int(cyclic_period)
        period >= 1 || throw(DomainError(period,
            "`cyclic_period` must be a positive integer."))
        n % period == 0 || throw(ArgumentError(
            "`cyclic_period` must divide the code length $n."))
        # Lifted circulants from F_q[x]/(x^ℓ-1) are invariant under the
        # simultaneous rotation of every ℓ-column block. Any nonzero word
        # has a representative with a 1 in the t = 0 slice.
        @constraint(model, sum(x[t] for t in 1:period:n) >= 1)
        verbose && println("Added cyclic slice cut of period $period.")
    end
    @objective(model, Min, sum(x))

    verbose && println("Added $parity_cut_count odd-set parity cuts.")
    threads > 0 && verbose && println("HiGHS threads = $threads.")
    optimize!(model)
    status = termination_status(model)
    if status == JuMP.MOI.INFEASIBLE
        return -1, zeros(Int, n), :infeasible
    elseif status == JuMP.MOI.TIME_LIMIT
        if has_values(model)
            witness = round.(Int, value.(x))
            return sum(witness), witness, :time_limit
        end
        return -1, zeros(Int, n), :time_limit
    elseif status != JuMP.MOI.OPTIMAL
        incumbent = has_values(model) ? objective_value(model) : missing
        bound = objective_bound(model)
        error("CSS minimum-distance ILP did not finish exactly " *
              "(status: $status, incumbent: $incumbent, bound: $bound).")
    end

    witness = round.(Int, value.(x))
    return sum(witness), witness, :optimal
end
