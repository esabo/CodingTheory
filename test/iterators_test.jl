@testitem "iterators.jl" begin 
    let
        len = 5
        weight = 2
        itr = CodingTheory.SubsetGrayCode(len, weight) 
        all_subsets_gray = Vector()
        v = collect(1: weight)
        rank = 1
        init_inds = [-1, -1, -1]
        state = (v, rank, init_inds)
        for i in 1:length(itr)
            push!(all_subsets_gray, deepcopy(v)) 
            next = iterate(itr, state)
            next === nothing && break
            (_, state) = next
            (v, _, _) = state
        end
    
        all_subsets = collect(CodingTheory.Oscar.subsets(collect(1:len), weight))
        @test length(all_subsets_gray) == length(all_subsets)
        @test all([u in all_subsets_gray for u in all_subsets])
        @test all([u in all_subsets for u in all_subsets_gray])
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