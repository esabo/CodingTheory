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

@time _modified_hawick_james(D, L.n)
