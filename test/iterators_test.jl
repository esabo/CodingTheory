@testitem "iterators.jl" begin 
    let
        v_len = 15
        v_weight = 7
        sgc = CodingTheory.SubsetGrayCode(v_len, v_weight) 
        all_subsets_gray = fill(fill(0, v_weight), length(sgc))
        vec = collect(1:sgc.k)
        state = (vec, 1, Array{Int}([-1, -1, -1]))
        for i in 1:length(sgc)
            all_subsets_gray[i] = deepcopy(vec)
                next = iterate(sgc, state)
                next === nothing && break
                (_, state) = next
                (vec, _, _) = state
        end
        sort!(all_subsets_gray)
        all_subsets_hecke = CodingTheory.Hecke.subsets(collect(1:v_len), v_weight)
        sort!(all_subsets_hecke)
        @assert all_subsets_gray == all_subsets_hecke 
    end

    function pushOrDel(a::Set{T}, b::T...) where T <: Any
        c = deepcopy(a)
        pushOrDel!(c,b)
        return c
    end

    function pushOrDel!(a::Set{T}, b::T...) where T <: Any
        for elem in b
            if elem in a
                delete!(a, elem)
            else
                push!(a, elem)
            end
        end
        return nothing
    end

    let 
        tuple0 = Vector{UInt}([1,2,3])
        kk = UInt(3)
        nn = UInt(5)
        rank = CodingTheory._kSubsetRevDoorRank(tuple0, kk)
        @test rank == 0
        result0 = zeros(UInt, kk)
        CodingTheory._kSubsetRevDoorUnrank(rank, nn, result0)
        @test result0 == tuple0
  
        tuple0 = UInt.([1,3,5])
        kk = UInt(3)
        nn = UInt(5)
        rank = CodingTheory._kSubsetRevDoorRank(tuple0, kk)
        @test rank == 7
        result1 = zeros(UInt, kk)
        CodingTheory._kSubsetRevDoorUnrank(rank, nn, result1)
        result1 = Int.(result1)
        @test result1 == [1, 3, 5] 
  
        kk = UInt(3)
        nn = UInt(5)
        result0 = zeros(UInt, kk)
        results = Set()
        bin = binomial(nn, kk) 
        for i::BigInt in collect(0: bin-1)
          CodingTheory._kSubsetRevDoorUnrank(i, nn, result0)
          pushOrDel!(results, deepcopy(result0))
        end
        @test length(results) == bin 
    end
end