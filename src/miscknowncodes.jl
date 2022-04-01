include("linearcode.jl")

function repetitioncode(n::Integer)
    F, _ = FiniteField(2, 1, "α")
    G = matrix(F, 1, n, [1 for i in 1:n])
    M2 = MatrixSpace(F, n - 1, 1)
    M3 = MatrixSpace(F, n - 1, n - 1)
    H = hcat(M2([1 for i in 1:(n - 1)]), M3(1))
    return LinearCode(F, n, 1, n, G, G, H, H, G, H, missing)
end

# function Hammingcode(p::Integer, r::Integer)
#
# end

# function tetracode()
#     return Hammingcode(3, 2)
# end

# both of these are actually the extended Golay codes
function Golay(p::Integer)
    if p == 2
        F, _ = FiniteField(2, 1, "α")
        M = MatrixSpace(F, 12 , 12)
        A = M([0 1 1 1 1 1 1 1 1 1 1 1;
             1 1 1 0 1 1 1 0 0 0 1 0;
             1 1 0 1 1 1 0 0 0 1 0 1;
             1 0 1 1 1 0 0 0 1 0 1 1;
             1 1 1 1 0 0 0 1 0 1 1 0;
             1 1 1 0 0 0 1 0 1 1 0 1;
             1 1 0 0 0 1 0 1 1 0 1 1;
             1 0 0 0 1 0 1 1 0 1 1 1;
             1 0 0 1 0 1 1 0 1 1 1 0;
             1 0 1 0 1 1 0 1 1 1 0 0;
             1 1 0 1 1 0 1 1 1 0 0 0;
             1 0 1 1 0 1 1 1 0 0 0 1])
        G = hcat(M(1), A)
        H = hcat(-A', M(1))
        return LinearCode(F, 24, 12, 8, G, G, H, H, G, H, missing)
    elseif p == 3
        F, _ = FiniteField(3, 1, "α")
        M = MatrixSpace(F, 6 , 6)
        A = M([0 1 1 1 1 1;
               1 0 1 -1 -1 1;
               1 1 0 1 -1 -1;
               1 -1 1 0 1 -1;
               1 -1 -1 1 0 1;
               1 1 -1 -1 1 0])
        G = hcat(M(1), A)
        H = hcat(-A', M(1))
        return LinearCode(F, 12, 6, 6, G, G, H, H, G, H, missing)
    else
        error("Golay code not implemented for q = $q.")
    end
end
