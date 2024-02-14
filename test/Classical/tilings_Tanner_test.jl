@testset "Classical/tilings.jl & Classical/Tanner.jl" begin
    using Graphs
    CodingTheory.Oscar.GAP.Packages.load("LINS");

    min_index = 250;
    max_index = 5000;
    # F, _ = finite_field(2)
    F = GF(2)

    # first test case
    local_code = HammingCode(2, 3);
    H = parity_check_matrix(local_code)
    locs_wts = Vector{Int}()
    for i in 1:nrows(H)
        push!(locs_wts, wt(H[i, :]))
    end

    g = r_s_group(3, 7);
    subgroups = normal_subgroups(g, max_index)
    for subgroup in subgroups
        # for this test case, the subgroup are numbers [3, 4, 5, 6, 7, 8, 9]
        if is_fixed_point_free(subgroup, g) && GAP.Globals.Index(g.group, subgroup) > min_index
            adj = sparse(transpose(coset_intersection([2, 3], [1, 3], subgroup, g)))
            code = Tanner_code(adj, local_code)
            # code = LinearCode(matrix(F, code), true) # remove later
            @test code.k >= adj.n - adj.m * (local_code.n - local_code.k)

            flag = true
            for i in 1:nrows(code.H)
                wt(code.H[i, :]) ∈ locs_wts || (flag = false;)
            end
            @test flag
        end
    end

    # second test case
    # C1 = GAP.Globals.BestKnownLinearCode(5, 2, GAP.Globals.GF(2))
    # x = GAP.Globals.GeneratorMat(C1)
    # y = [GAP.Globals.Int(x[i, j]) for i in 1:2, j in 1:5]
    y = [0 0 1 1 1; 1 1 0 1 1]
    z = matrix(F, y)
    C_loc = LinearCode(z)
    G_test = Graphs.complete_graph(6)
    EVI = sparse(transpose(Graphs.incidence_matrix(G_test)))
    H1 = Tanner_code(EVI, C_loc)
    EVI_G, left, right = edge_vertex_incidence_graph(G_test)
    H2 = Tanner_code(EVI_G, left, right, C_loc)
    @test parity_check_matrix(H1) == parity_check_matrix(H2)

    # third test case
    g = r_s_group(5, 4)
    test_mat = CodingTheory.CoxeterMatrix(3, [1, 5, 2, 1, 4, 1])
    @test g.cox_mat == test_mat
    subgroups = normal_subgroups(g, 100)
    @test length(subgroups) == 7
    filter!(x -> is_fixed_point_free(x, g), subgroups)
    @test length(subgroups) == 0

    # fourth test case
    g = r_s_group(8, 3)
    test_mat = CodingTheory.CoxeterMatrix(3, [1, 8, 2, 1, 3, 1])
    @test g.cox_mat == test_mat
    subgroups = normal_subgroups(g, 100)
    @test length(subgroups) == 10
    filter!(x -> is_fixed_point_free(x, g), subgroups)
    @test length(subgroups) == 1

    # fifth test case
    g = triangle_group(4, 3, 3)
    test_mat = CodingTheory.CoxeterMatrix(3, [1, 4, 3, 1, 3, 1])
    @test g.cox_mat == test_mat
    subgroups = normal_subgroups(g, 100)
    @test length(subgroups) == 6
    filter!(x -> is_fixed_point_free(x, g), subgroups)
    @test length(subgroups) == 2

    # sixth test case
    g = q_r_s_group(4, 3, 5)
    test_mat = CodingTheory.CoxeterMatrix(4, [1, 4, 2, 2, 1, 3, 2, 1, 5, 1])
    @test g.cox_mat == test_mat
    subgroups = normal_subgroups(g, 100)
    @test length(subgroups) == 6
    filter!(x -> is_fixed_point_free(x, g), subgroups)
    @test length(subgroups) == 0

    # seventh test case
    g = star_tetrahedron_group(5, 3, 3)
    test_mat = CodingTheory.CoxeterMatrix(4, [1, 5, 3, 3, 1, 2, 2, 1, 2, 1])
    @test g.cox_mat == test_mat
    subgroups = normal_subgroups(g, 100)
    @test length(subgroups) == 3
    filter!(x -> is_fixed_point_free(x, g), subgroups)
    @test length(subgroups) == 0
end
