include("linearcode.jl")

function repetitioncode(n::Integer)
    M = MatrixSpace(Nemo.GF(2), 1, n)
    G = M([1 for i in 1:n])
    M2 = MatrixSpace(Nemo.GF(2), n - 1, 1)
    M3 = MatrixSpace(Nemo.GF(2), n - 1, n - 1)
    H = hcat(M2([1 for i in 1:(n - 1)]), M3(1))
    return LinearCode(2, n, 1, G, G, H, G, H)
end

function Hammingcode(p::Integer, r::Integer)

end

function tetracode()
    return Hammingcode(3, 2)
end

function Golay(p::Integer)
    if p == 2
        M = MatrixSpace(Nemo.GF(2), 12 , 12)
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
        return LinearCode(2, 24, 12, G, G, H, G, H)
    elseif p == 3
        M = MatrixSpace(Nemo.GF(3), 6 , 6)
        A = M([0 1 1 1 1 1;
               1 0 1 -1 -1 1;
               1 1 0 1 -1 -1;
               1 -1 1 0 1 -1;
               1 -1 -1 1 0 1;
               1 1 -1 -1 1 0])
        G = hcat(M(1), A)
        H = hcat(-A', M(1))
        return LinearCode(3, 12, 6, G, G, H, G, H)
    else
        error("Golay code not implemented for q = $q.")
    end
end
