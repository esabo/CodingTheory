"""
    tannercode(vertexedgeadj::SparseMatrixCSC{Int,Int}, localcode::LinearCode)

Return the tanner code constructed from the graph defined by `vertexedgeadj` and `localcode`.
"""
function tannercode(vertexedgeadj::SparseMatrixCSC{Int,Int}, localcode::LinearCode)
    localcode.n == length(vertexedgeadj[1,:].nzind) || throw(ArgumentError("Local code length must equal graph valency."))
    n = vertexedgeadj.n
    m = size(vertexedgeadj)[1] * (localcode.n - localcode.k)
    # println(m, ",", n)
    F, _ = FiniteField(2, 1, "α")
    pcm = zero_matrix(F, m, n)
    l = 1
    for i in 1:size(vertexedgeadj)[1]
        edges = vertexedgeadj[i,:].nzind
        for c in 1:size(localcode.H)[1]
            check = vec(localcode.H[c,:] .!= F(0))
            # println(check)
            supp = getindex(edges, check)
            # println(supp)
            for idx in supp
                pcm[l,idx] = F(1)
            end
            l += 1
        end
    end
    LinearCode(pcm, true)
end
