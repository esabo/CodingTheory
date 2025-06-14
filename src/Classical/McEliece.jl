


function McEliece_attack(
    G::CTMatrix,
    w::CTMatrix,
    t::Int;
    max_iters::Int = 10000,
    verbose::Bool = false,
)

    # check 0 matrix, 0 dimensions

    # G - gen matrix
    # w - received vector
    # t - # errors in w
    # max_iters - maximum number of iterations
    n = ncols(G)
    G2 = _Flint_matrix_to_Julia_T_matrix(G2, UInt8)
    w2 = _Flint_matrix_to_Julia_T_matrix(w, UInt8)

    num_thrds = Threads.nthreads()
    verbose && println("Detected $num_thrds threads.")
    thread_load = Int(floor(max_iters / num_thrds))
    remaining = max_iters - thread_load * num_thrds
    flag = Threads.Atomic{Bool}(true)
    answers = [zeros(Int, 1, n) for _ = 1:num_thrds]
    which_ans = [false for _ = 1:num_thrds]
    perms = [collect(1:n) for _ = 1:num_thrds]

    verbose && (prog_meter = Progress(max_iters);)
    Threads.@threads for th = 1:num_thrds
        G_loc = similar(G2)
        σ_loc = shuffle(1:n)
        _col_permutation!(G2, G_loc, σ_loc)
        # doesn't work for ints yet
        # type stability with missing
        # maybe directly return perm in array style
        k, P_loc = _rref_col_swap!(G_loc)
        G_loc = _remove_empty(G_loc, :rows)
        # redo multiplication
        σ_loc *= P_loc
        w_loc = similar(w2)
        w_loc = _col_permutation!(w2, w_loc, σ_loc)

        for _ = 1:(thread_load+(th<=remaining ? 1 : 0))
            if flag[]
                wt = 0
                @inbounds for j = (k+1):n
                    isodd(w_loc[i, j]) && (wt += 1;)
                end
                if wt == t
                    which_ans[th] = true
                    answers[th][1, (k+1):n] .= w_loc[1, (k+1):n]
                    answers[th][1, invperm(σ_loc)]
                    perms[th] .= σ_loc
                    Threads.atomic_cas!(flag, true, false)
                    break
                else
                    shuffle!(σ_loc)
                    _col_permutation!(G2, G_loc, σ_loc)
                    _, P_loc = _rref_col_swap!(G_loc)
                    w_loc = _col_permutation!(w2, w_loc, σ_loc)
                    σ_loc *= P_loc
                end
            else
                break
            end
            verbose && next!(prog_meter)
        end
    end
    verbose && finish!(prog_meter)

    i = findfirst(which_ans, true)
    # cheaper to not return the matrix here
    return matrix(base_ring(G), answers[i]), perms[i]
end

function Lee_Brickell_attack(
    G::CTMatrix,
    w::CTMatrix,
    t::Int,
    p::Int;
    max_iters::Int = 10000,
    verbose::Bool = false,
)
    # check 0 matrix, 0 dimensions

    # G - gen matrix
    # w - received vector
    # t - # errors in w
    # max_iters - maximum number of iterations
    n = ncols(G)
    G2 = _Flint_matrix_to_Julia_T_matrix(G2, UInt8)
    w2 = _Flint_matrix_to_Julia_T_matrix(w, UInt8)

    num_thrds = Threads.nthreads()
    verbose && println("Detected $num_thrds threads.")
    thread_load = Int(floor(max_iters / num_thrds))
    remaining = max_iters - thread_load * num_thrds
    flag = Threads.Atomic{Bool}(true)
    answers = [zeros(Int, 1, n) for _ = 1:num_thrds]
    which_ans = [false for _ = 1:num_thrds]
    perms = [collect(1:n) for _ = 1:num_thrds]

    # update size for p and size of GrayCode
    verbose && (prog_meter = Progress(max_iters);)
    Threads.@threads for th = 1:num_thrds
        G_loc = similar(G2)
        σ_loc = shuffle(1:n)
        _col_permutation!(G2, G_loc, σ_loc)
        # doesn't work for ints yet
        # type stability with missing
        # maybe directly return perm in array style
        k, P_loc = _rref_col_swap!(G_loc)
        G_loc = _remove_empty(G_loc, :rows)
        # redo multiplication
        σ_loc *= P_loc
        w_loc = similar(w2)
        w_loc = _col_permutation!(w2, w_loc, σ_loc)
        q_A = zeros(Int, 1, n - k)
        for _ = 1:(thread_load+(th<=remaining ? 1 : 0))
            if flag[]
                wt = 0
                @inbounds for j = (k+1):n
                    isodd(w_loc[i, j]) && (wt += 1;)
                end

                if wt ≤ t
                    which_ans[th] = true
                    answers[th][1, (k+1):n] .= w_loc[1, (k+1):n]
                    answers[th][1, invperm(σ_loc)]
                    perms[th] .= σ_loc
                    Threads.atomic_cas!(flag, true, false)
                    break
                end

                for r = 1:p
                    if flag[]
                        for m in GrayCode(k, r)
                            if flag[]
                                LinearAlgebra.mul!(q_A, m, view(G_loc, :, (k+1):n))
                                q_A .+= w_loc[1, (k+1):n]
                                wt = 0
                                @inbounds for j in axes(q_A, 2)
                                    isodd(q_A[1, j]) && (wt += 1;)
                                end

                                if wt ≤ t - r
                                    which_ans[th] = true
                                    answers[th][1, 1:k] .= m
                                    answers[th][1, (k+1):n] .= q_A
                                    answers[th][1, invperm(σ_loc)]
                                    perms[th] .= σ_loc
                                    Threads.atomic_cas!(flag, true, false)
                                    break
                                end
                            else
                                break
                            end
                        end
                    else
                        break
                    end
                end

                if flag[]
                    shuffle!(σ_loc)
                    _col_permutation!(G2, G_loc, σ_loc)
                    _, P_loc = _rref_col_swap!(G_loc)
                    w_loc = _col_permutation!(w2, w_loc, σ_loc)
                    σ_loc *= P_loc
                    verbose && next!(prog_meter)
                else
                    break
                end
            else
                break
            end
        end
    end
    verbose && finish!(prog_meter)

    i = findfirst(which_ans, true)
    # cheaper to not return the matrix here
    return matrix(base_ring(G), answers[i]), perms[i]
end

function Leon_attack(
    G::CTMatrix,
    w::CTMatrix,
    t::Int,
    p::Int,
    l::Int;
    max_iters::Int = 10000,
    verbose::Bool = false,
)

    # G - gen matrix
    # w - received vector
    # t - # errors in w
    # max_iters - maximum number of iterations
    n = ncols(G)
    G2 = _Flint_matrix_to_Julia_T_matrix(G2, UInt8)
    w2 = _Flint_matrix_to_Julia_T_matrix(w, UInt8)

    num_thrds = Threads.nthreads()
    verbose && println("Detected $num_thrds threads.")
    thread_load = Int(floor(max_iters / num_thrds))
    remaining = max_iters - thread_load * num_thrds
    flag = Threads.Atomic{Bool}(true)
    answers = [zeros(Int, 1, n) for _ = 1:num_thrds]
    which_ans = [false for _ = 1:num_thrds]
    perms = [collect(1:n) for _ = 1:num_thrds]

    # update size for p and size of GrayCode
    verbose && (prog_meter = Progress(max_iters);)
    Threads.@threads for th = 1:num_thrds
        G_loc = similar(G2)
        σ_loc = shuffle(1:n)
        _col_permutation!(G2, G_loc, σ_loc)
        # doesn't work for ints yet
        # type stability with missing
        # maybe directly return perm in array style
        k, P_loc = _rref_col_swap!(G_loc)
        G_loc = _remove_empty(G_loc, :rows)
        # redo multiplication
        σ_loc *= P_loc
        w_loc = similar(w2)
        w_loc = _col_permutation!(w2, w_loc, σ_loc)
        q_P = zeros(Int, 1, l)
        q_B = zeros(Int, 1, n - k - l)
        for _ = 1:(thread_load+(th<=remaining ? 1 : 0))
            if flag[]
                wt = 0
                @inbounds for j = (k+1):n
                    isodd(w_loc[i, j]) && (wt += 1;)
                end

                if wt ≤ t
                    which_ans[th] = true
                    answers[th][1, (k+1):n] .= w_loc[1, (k+1):n]
                    answers[th][1, invperm(σ_loc)]
                    perms[th] .= σ_loc
                    Threads.atomic_cas!(flag, true, false)
                    break
                end

                for r = 1:p
                    if flag[]
                        for m in GrayCode(k, r)
                            if flag[]
                                LinearAlgebra.mul!(q_P, m, view(G_loc, :, (k+1):(k+l)))
                                q_P .+= w_loc[1, (k+1):(k+l)]
                                wt_P = 0
                                @inbounds for j in axes(q_P, 2)
                                    isodd(q_P[1, j]) && (wt_P += 1;)
                                end

                                if wt_P ≤ p - r
                                    LinearAlgebra.mul!(q_B, m, view(G_loc, :, (k+l+1):n))
                                    q_B .+= w_loc[1, (k+l+1):n]
                                    wt_B = 0
                                    @inbounds for j in axes(q_B, 2)
                                        isodd(q_B[1, j]) && (wt_B += 1;)
                                    end

                                    if wt_B ≤ t - r - wt_p
                                        which_ans[th] = true
                                        answers[th][1, 1:k] .= m
                                        answers[th][1, (k+1):(k+l)] .= q_P
                                        answers[th][1, (k+l+1):n] .= q_B
                                        answers[th][1, invperm(σ_loc)]
                                        perms[th] .= σ_loc
                                        Threads.atomic_cas!(flag, true, false)
                                        break
                                    end
                                end
                            else
                                break
                            end
                        end
                    else
                        break
                    end
                end

                if flag[]
                    shuffle!(σ_loc)
                    _col_permutation!(G2, G_loc, σ_loc)
                    _, P_loc = _rref_col_swap!(G_loc)
                    w_loc = _col_permutation!(w2, w_loc, σ_loc)
                    σ_loc *= P_loc
                    verbose && next!(prog_meter)
                else
                    break
                end
            else
                break
            end
        end
    end
    verbose && finish!(prog_meter)

    i = findfirst(which_ans, true)
    # cheaper to not return the matrix here
    return matrix(base_ring(G), answers[i]), perms[i]
end


# my other is called Sterns_attack so won't mess it up and this file isn't imported yet
function generalized_Stern_attack(
    G::CTMatrix,
    w::CTMatrix,
    t::Int,
    p::Int,
    l::Int,
    s::Int;
    max_iters::Int = 10000,
    verbose::Bool = false,
)

    n = ncols(G)
    G2 = _Flint_matrix_to_Julia_T_matrix(G2, UInt8)
    w2 = _Flint_matrix_to_Julia_T_matrix(w, UInt8)

    num_thrds = Threads.nthreads()
    verbose && println("Detected $num_thrds threads.")
    thread_load = Int(floor(max_iters / num_thrds))
    remaining = max_iters - thread_load * num_thrds
    flag = Threads.Atomic{Bool}(true)
    answers = [zeros(Int, 1, n) for _ = 1:num_thrds]
    which_ans = [false for _ = 1:num_thrds]
    perms = [collect(1:n) for _ = 1:num_thrds]

    # update size for p and size of GrayCode
    verbose && (prog_meter = Progress(max_iters);)
    Threads.@threads for th = 1:num_thrds
        G_loc = similar(G2)
        σ_loc = shuffle(1:n)
        _col_permutation!(G2, G_loc, σ_loc)
        # doesn't work for ints yet
        # type stability with missing
        # maybe directly return perm in array style
        k, P_loc = _rref_col_swap!(G_loc)
        G_loc = _remove_empty(G_loc, :rows)
        # redo multiplication
        σ_loc *= P_loc
        w_loc = similar(w2)
        w_loc = _col_permutation!(w2, w_loc, σ_loc)
        q_i = zeros(Int, 1, l)
        q_sum = zeros(Int, 1, l)
        q_B = zeros(Int, 1, n - k - l)
        m_B = zeros(Int, 1, k)

        nr_Xi = Int(floor(k / s))
        nr_Y = k - s * nr_Xi

        # first check none of these submatrices are 0

        Xis = [Dict{Matrix{UInt8},Vector{Matrix{UInt8}}}() for _ = 1:(s-1)]
        size_guess = binomial(BigInt(nr_Xi), p)
        for i = 1:(s-1)
            sizehint!(Xis[i], size_guess)
        end

        for _ = 1:(thread_load+(th<=remaining ? 1 : 0))
            if flag[]
                # does this check and if statement still hold here?
                wt = 0
                @inbounds for j = (k+1):n
                    isodd(w_loc[i, j]) && (wt += 1;)
                end

                if wt ≤ t
                    which_ans[th] = true
                    answers[th][1, (k+1):n] .= w_loc[1, (k+1):n]
                    answers[th][1, invperm(σ_loc)]
                    perms[th] .= σ_loc
                    Threads.atomic_cas!(flag, true, false)
                    break
                end

                for m in GrayCode(nr_Xi, p)
                    for i = 1:(s-1)
                        LinearAlgebra.mul!(
                            q_i,
                            m,
                            view(G_loc, ((i-1)*nr_Xi+1):(i*nr_Xi), (k+1):(k+l)),
                        )
                        q_i .+= w_loc[1, (k+1):(k+l)]
                        if haskey(Xis[i], q_i)
                            push(Xis[i][q_i], m)
                        else
                            Xis[i][q_i] = m
                        end
                    end
                end

                # allocates a pointer, don't ask for dynamically in a loop
                keys_Xis = [keys(Xi_s) for _ = 1:length(Xis)]
                vals_Xis = [values(Xi_s) for _ = 1:length(Xis)]
                for m in GrayCode(nr_Y, p)
                    if flag[]
                        # q_i now representing q_Y
                        LinearAlgebra.mul!(q_i, m, view(G_loc, (k-nr_Y+1):k, (k+1):(k+l)))
                        q_i .+= w_loc[1, (k+1):(k+l)]
                        for iter_q in Nemo.AbstractAlgebra.ProductIterator(
                            [1:length(Xis) for _ = 1:(s-1)],
                            inplace = true,
                        )
                            if flag[]
                                q_sum .= keys_Xis[1][iter_q[1]]
                                for i = 2:length(iter_q)
                                    q_sum .+= keys_Xis[i][iter_q[i]]
                                end

                                if q_sum == q_i
                                    for iter_m in Nemo.AbstractAlgebra.ProductIterator([
                                        1:length(
                                            vals_Xis[i][iter_q[i]] for i = 1:length(Xis)
                                        ),
                                    ])
                                        # make m_B, q_B
                                        for i = 1:(s-1)
                                            m_B[1, ((i-i)*nr_Xi+1):(i*nr_Xi)] .=
                                                vals_Xis[i][iter[i]][iter_m[i]]
                                        end
                                        LinearAlgebra.mul!(
                                            q_B,
                                            m_B,
                                            view(G_loc, :, (k+l+1):n),
                                        )
                                        q_B .+= w_loc[1, (k+l+1):n]
                                        wt_B = 0
                                        @inbounds for j in axes(q_B, 2)
                                            isodd(q_B[1, j]) && (wt_B += 1;)
                                        end

                                        if wt_B ≤ t - 2p
                                            which_ans[th] = true
                                            answers[th][1, 1:k] .= m_B
                                            answers[th][1, (k+l+1):n] .= q_B
                                            answers[th][1, invperm(σ_loc)]
                                            perms[th] .= σ_loc
                                            Threads.atomic_cas!(flag, true, false)
                                            break
                                        end
                                    end
                                end
                            else
                                break
                            end
                        end
                    else
                        break
                    end
                end

                if flag[]
                    shuffle!(σ_loc)
                    _col_permutation!(G2, G_loc, σ_loc)
                    _, P_loc = _rref_col_swap!(G_loc)
                    w_loc = _col_permutation!(w2, w_loc, σ_loc)
                    σ_loc *= P_loc
                    verbose && next!(prog_meter)
                else
                    break
                end
            else
                break
            end
        end
    end
    verbose && finish!(prog_meter)

    i = findfirst(which_ans, true)
    # cheaper to not return the matrix here
    return matrix(base_ring(G), answers[i]), perms[i]
end
Stern_attack(
    G::CTMatrix,
    w::CTMatrix,
    t::Int,
    p::Int,
    l::Int,
    s::Int;
    max_iters::Int = 10000,
    verbose::Bool = false,
) = generalized_Stern_attack(G, w, t, p, l, 1, max_iters = max_iters, verbose = verbose)

function Canteaut_Chabaud_attack()

end
