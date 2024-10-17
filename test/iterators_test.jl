@testitem "iterators.jl" begin 
  @testset "Subset Gray Code" begin
    len1=5
    weight1=2
    gi=CodingTheory.SubsetGrayCode(len1, weight1) # gray iterator 
    i = 0
    all_subsets_gray = Vector()
    my_vec = collect(1: weight1)
    my_state = (my_vec, 1, Array{Int}([-1,-1,-1]))
    for i in 1:length(gi)
        push!(all_subsets_gray, deepcopy(my_vec)) 
        next = iterate(gi, my_state)
        next === nothing && break
        (_, my_state) = next
        (my_vec, _, _) = my_state
    end

    subsets = CodingTheory.Oscar.subsets(collect(1:len1), weight1)
    all_subsets = collect(subsets)
    @test length(all_subsets_gray) == length(all_subsets)
    @test all([u in all_subsets_gray for u in all_subsets])
    @test all([u in all_subsets for u in all_subsets_gray])
  end
end
