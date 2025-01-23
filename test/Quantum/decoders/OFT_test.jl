@testitem "Quantum/decoders/OTF.jl" begin
    using Random, Oscar, CodingTheory

    @testset "Ordered Tanner Forest" begin
        # first define Ton's Julia transpilation of the original Python
        function get_indices(H::Matrix{UInt8}, ordered_indices::Vector{Int} = collect(1:size(H, 2)))
            # We generate the linked_list
            m = size(H, 1)
            linked_list = ones(Int, 2 * m)
            for i in 1:m
                linked_list[2 * i - 1] = i
            end
        
            # Initialize an to hold the results, the True elements will have their probabilities updated to a very low value.
            result_indices = falses(size(H, 2))
        
            # Iterate over the columns in the specified order
            for col in ordered_indices
                # println("col $col")
                # Get the rows that have True in the current column
                # TODO, this can be simplified by using a sparse matrix instead of a PCM, following line is too expensive.
                true_rows = findall(row -> isone(H[row, col]), 1:size(H, 1))
                
                adding_index = false
                growing_depth = false
                row_list = Int[]
                max_depth = 0
                root = -1
        
                for true_row in true_rows
                    index = find(linked_list, true_row)
                    if !(index in row_list)
                        push!(row_list, index)
                        if linked_list[2 * index] > max_depth
                            root = index
                            growing_depth = false
                            max_depth = linked_list[2 * index]
                        elseif linked_list[2 * index] == max_depth
                            growing_depth = true
                        end
                    else
                        adding_index = true
                        break
                    end
                end
        
                if adding_index
                    result_indices[col] = true
                else
                    for index in row_list
                        linked_list[2 * index - 1] = root
                    end
                end
                if growing_depth
                    linked_list[2 * root] += 1
                end
                # println(linked_list)
            end
        
            return result_indices
        end

        function find(linked_list::Vector{Int}, index::Int)
            while linked_list[2 * index - 1] != index
                index = linked_list[2 * index - 1]
            end
            return index
        end

        # check it against the code in this library
        for S in [SteaneCode(), Q15RM(), RotatedSurfaceCode(3), RotatedSurfaceCode(5), RotatedSurfaceCode(7), GrossCode()]
            H_Int = CodingTheory._Flint_matrix_to_Julia_T_matrix(X_stabilizers(S), UInt8);
            ordering = shuffle(1:S.n);
            var_adj_list = [Int[] for _ in 1:size(H_Int, 2)];
            for r in 1:size(H_Int, 1)
                for c in 1:size(H_Int, 2)
                    if !iszero(H_Int[r, c])
                        push!(var_adj_list[c], r)
                    end
                end
            end
            E = sort(CodingTheory._select_erased_columns(H_Int, ordering, var_adj_list));
            T = findall(x -> x == true, get_indices(H_Int, ordering));
            @test E == T
            # test = E == T
            # println(test)
            # if !test
            #     println(ordering)
            #     println(E)
            #     println(T)
            # end
        end
    end
end
