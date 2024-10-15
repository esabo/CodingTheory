function test_enumerates_all_subsets()
  glen=5
  weight=2
  gi=CodingTheory.SubsetGrayCode(glen, weight) # gray iterator 
  i = 0
  all_subsets_gray = Vector()
  vec = collect(1:gi.k)
  state = (vec, 1, Array{Int}([-1,-1,-1]))
  for i in 1:length(gi)
      push!(all_subsets_gray, deepcopy(vec)) 
      next = iterate(gi, state)
      next === nothing && break
      (_, state) = next
      (vec, _, _) = state
  end

  all_subsets = collect(subsets(collect(1:glen), weight))
  @test length(all_subsets_gray) == length(all_subsets)
  @test all([u in all_subsets_gray for u in all_subsets])
  @test all([u in all_subsets for u in all_subsets_gray])
end