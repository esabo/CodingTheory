@testitem "iterators.jl" begin 
    @testset "SubsetGrayCode iterates over all weight k subsets of {1,..,n}" begin
        len = 15
        weight = 7
        sgc = CodingTheory.SubsetGrayCode(len, weight) 
        all_subsets_gray = fill(fill(0, weight), length(sgc))
        subset = collect(1:sgc.k)
        state = (subset, 1, fill(-1, 3))
        for i in 1:length(sgc)
            all_subsets_gray[i] = deepcopy(subset)
                next = iterate(sgc, state)
                next === nothing && break
                (_, state) = next
                (subset, _, _) = state
        end
        sort!(all_subsets_gray)
        all_subsets_hecke = CodingTheory.Oscar.subsets(collect(1:len), weight)
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

    @testset "Rank/Unrank functions" begin
        subset1 = Vector{UInt}([1, 2, 3]) # a subset of weight 3
        k1 = UInt(3)
        n1 = UInt(5)
        rank1 = CodingTheory._subset_rank(subset1, k1)
        @test rank1 == 0
        result1 = zeros(UInt, k1)
        CodingTheory._subset_unrank!(rank1, n1, result1)
        @test result1 == subset1
  
        subset2 = UInt.([1, 3, 5])
        k2 = UInt(3)
        n2 = UInt(5)
        rank2 = CodingTheory._subset_rank(subset2, k2)
        @test rank2 == 7
        result2 = zeros(UInt, k2)
        CodingTheory._subset_unrank!(rank2, n2, result2)
        @test result2 == subset2 
  
        k3 = UInt(3)
        n3 = UInt(5)
        result3 = zeros(UInt, k3)
        results = Set()
        bin = binomial(n3, k3) 
        for i::BigInt in collect(0: bin - 1)
          CodingTheory._subset_unrank!(i, n3, result3)
          pushOrDel!(results, deepcopy(result3))
        end
        all_subsets_hecke = Set(CodingTheory.Hecke.subsets(collect(1:n3), Int(k3)))
        @test results == all_subsets_hecke 
    end
end