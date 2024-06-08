function _modified_hawick_james(g::DiGraph{Int}, num_var_nodes::Int)
    nvg = Graphs.nv(g)
    local_cycles = [Vector{Vector{Int}}() for _ in 1:num_var_nodes]
    Threads.@threads for i in 1:num_var_nodes
        B = [Vector{Int}() for _ in Graphs.vertices(g)]
        blocked = zeros(Bool, nvg)
        stack = Vector{Int}()
        keys_Dict = Dict{Vector{Int}, Bool}()
        _circuit_recursive!(g, Graphs.vertices(g)[i], Graphs.vertices(g)[i], blocked, B, stack, local_cycles[i], keys_Dict)
    end
    return reduce(vcat, local_cycles)
end

function _circuit_recursive!(g::DiGraph{Int}, v1::Int, v2::Int, blocked::Vector{Bool},
    B::Vector{Vector{Int}}, stack::Vector{Int}, cycles::Vector{Vector{Int}},
    keys_Dict::Dict{Vector{Int}, Bool})

    if length(stack) + 1 <= 10
        flag = false
        push!(stack, v2)
        blocked[v2] = true

        # just put this entire thing in the if statement
        Av = Graphs.outneighbors(g, v2)
        for w in Av
            (w < v1) && continue
            if w == v1
                if length(stack) > 2
                    cycle = copy(stack)
                    # println("cycle: $stack")
                    # println("also checking: ", [cycle[1]; reverse(cycle[2:end])])
                    if !haskey(keys_Dict, cycle) && !haskey(keys_Dict, [cycle[1]; reverse(
                        cycle[2:end])])

                        push!(cycles, cycle)
                        keys_Dict[cycle] = true
                        # println("is new")
                    # else
                    #     println("is old")
                    end
                end
                flag = true
            elseif !blocked[w]
                # bit-wise or
                flag |= _circuit_recursive!(g, v1, w, blocked, B, stack, cycles, keys_Dict)
            end
        end

        if flag
            _unblock!(v2, blocked, B)
        else
            for w in Av
                (w < v1) && continue
                if !(v2 in B[w])
                    push!(B[w], v2)
                end
            end
        end
        
        pop!(stack)
    else
        flag = true
    end
    return flag
end

function _unblock!(v::Int, blocked::Vector{Bool}, B::Vector{Vector{Int}})
    blocked[v] = false
    wPos = 1
    Bv = B[v]
    while wPos <= length(Bv)
        w = Bv[wPos]
        old_length = length(Bv)
        filter!(v -> v != w, Bv)
        wPos += 1 - (old_length - length(Bv))
        if blocked[w]
            _unblock!(w, blocked, B)
        end
    end
    return nothing
end

# @time _modified_hawick_james(D, L.n)

function remove_cycles(H::Union{Matrix{<: Integer}, CTMatrixTypes}, n_max::Int)
    n_max ≥ 4 || throw(DomainError("n_max must be an even number greater than four"))
    iseven(n_max) || (n_max -= 1;)
    size(H, 1) ≥ 1 && size(H, 2) ≥ 1 || throw(ArgumentError("Matrix cannot have a zero dimension"))

    # also used in Tanner.jl, perhaps should make into its own function
    # A = _adjacency_matrix_from_code(H)
    typeof(H) <: CTMatrixTypes ? (I = _Flint_matrix_to_Julia_int_matrix(H);) : (I = H;)
    nr_H, nc_H = size(I)
    A = vcat(hcat(zeros(Int, nc_H, nc_H), transpose(I)), hcat(I, zeros(Int, nr_H, nr_H)))
    # display(A)
    matrices = [A^i for i in 0:n_max - 1]
    # for M in matrices
    #     display(M)
    #     println(" ")
    # end
    nr, nc = size(A)
    for n in 4:2:n_max
        n2 = div(n, 2)
        n2m1 = n2 - 1
        n2m2 = n2 - 2
        nm1 = n - 1
        for i in 1:nr
            for j in 1:nc
                if i ≠ j
                    if matrices[n2 + 1][i, j] ≥ 2 && iszero(matrices[n2m2 + 1][i, j])
                        # println("$i, $j")
                        # println(matrices[n2 + 1][i, j])
                        temp = Int[]
                        for k in 1:nc
                            # should be exactly two of these
                            if matrices[n2m1 + 1][i, k] > 0 && matrices[1 + 1][j, k] == 1
                                append!(temp, k)
                            end
                        end
                        
                        if !isempty(temp)
                            k = rand(temp)
                            C_e = Int[]
                            for x in 1:nc
                                if iszero(matrices[nm1 + 1][j, x]) && iszero(matrices[nm1 + 1][k, x])
                                    append!(C_e, x)
                                end
                            end

                            E_e = Vector{Tuple{Int, Int}}()
                            for v1 in C_e
                                for v2 in C_e
                                    if matrices[1 + 1][v1, v2] > 0
                                        push!(E_e, (v1, v2))
                                    end
                                end
                            end

                            # this handles the case that C_e is empty
                            if !isempty(E_e)
                                (l, m) = rand(E_e)
                                # println("here")
                                # remove two old edges
                                matrices[1 + 1][j, k] -= 1
                                matrices[1 + 1][k, j] -= 1
                                matrices[1 + 1][l, m] -= 1
                                matrices[1 + 1][m, l] -= 1

                                # add two new edges
                                matrices[1 + 1][j, m] += 1
                                matrices[1 + 1][m, j] += 1
                                matrices[1 + 1][l, k] += 1
                                matrices[1 + 1][k, l] += 1

                                # update matrices in indices j, k, l, m
                                for indx in 2:n_max - 1
                                    for c in 1:nc
                                        matrices[indx + 1][j, c] = dot(view(matrices[indx - 1 + 1], j, :), view(matrices[1 + 1], :, c))
                                        matrices[indx + 1][k, c] = dot(view(matrices[indx - 1 + 1], k, :), view(matrices[1 + 1], :, c))
                                        matrices[indx + 1][l, c] = dot(view(matrices[indx - 1 + 1], l, :), view(matrices[1 + 1], :, c))
                                        matrices[indx + 1][m, c] = dot(view(matrices[indx - 1 + 1], m, :), view(matrices[1 + 1], :, c))
                                    end

                                    for r in 1:nc
                                        matrices[indx + 1][r, j] = dot(view(matrices[indx - 1 + 1], r, :), view(matrices[1 + 1], :, j))
                                        matrices[indx + 1][r, k] = dot(view(matrices[indx - 1 + 1], r, :), view(matrices[1 + 1], :, k))
                                        matrices[indx + 1][r, l] = dot(view(matrices[indx - 1 + 1], r, :), view(matrices[1 + 1], :, l))
                                        matrices[indx + 1][r, m] = dot(view(matrices[indx - 1 + 1], r, :), view(matrices[1 + 1], :, m))
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return matrices[1 + 1][nc_H + 1:end, 1:nc_H]
end
# remove_cycles(L::LDPCCode, n_max::Int) = LDPCCode(remove_cycles(parity_check_matrix(L), n_max))
remove_cycles(L::LDPCCode, n_max::Int) = remove_cycles(parity_check_matrix(L), n_max)
