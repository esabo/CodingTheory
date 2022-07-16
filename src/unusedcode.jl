# const DATAFILE = joinpath(@__DIR__, "..", "deps", "conwaypolynomials.data")
const DATAFILE = "conwaypolynomials.data"
_conwaypolynomials = nothing

function wgetconwaypolynomials()
    url = "https://gist.githubusercontent.com/tkluck/e1cd1746c69aa17e4a37114d22649627/raw/7fbe9763fae27f14924262ad03606f1c3af4400e/CPImport.txt"
    ResDataType = Dict{Tuple{Int,Int}, Vector{Int}}
    result = ResDataType()

    data = String(HTTP.request("GET", url).body)
    data = replace(data, r"allConwayPolynomials := " => "")
    data = replace(data, r"];"                       => "]")

    list = JSON.parse(data)
    # Frank Lübeck's data format signifies the end of the array with a null marker
    pop!(list)

    result = ResDataType((p,n) => coeffs for (p,n,coeffs) in list)

    open("conwaypolynomials.data", "w") do f
        serialize(f, result)
    end
end

# initialize conway table and then later but this into BCH constructor
function conwaypolynomial(p::Integer, n::Integer)
    # strip this top part out into above
    global _conwaypolynomials
    if _conwaypolynomials == nothing
        # try to do this, if error then do wgetconwaypolynomials and redo
        _conwaypolynomials = deserialize(open(DATAFILE))
    end
    if (p, n) ∉ keys(_conwaypolynomials)
        error("No Conway polynomial for $p^$n in the database.")
    end
    return _conwaypolynomials[p, n]
end

# slower than built-in command
function ⊗(A::T, B::T) where T <: Union{fq_nmod_mat, gfp_mat}
    if base_ring(A) != base_ring(B)
        error("Matrices in tensor product must be over the same ring.")
    end

    M = MatrixSpace(base_ring(A), size(A, 1) * size(B, 1), size(A, 2) * size(B, 2))
    C = M(0)
    if iszero(A) || iszero(B)
        return M(0)
    else
        (rB, cB) = size(B)
        for ar in 1:size(A, 1)
            for ac in 1:size(A, 2)
                C[(1 + (ar - 1) * rB):(ar * rB), (1 + (ac - 1) * cB):(ac * cB)] = A[ar, ac] * B
            end
        end
    end

    return C
end

function Base.:(+)(a::Vector{Char}, b::Vector{Char})
    if length(a) != length(b)
        error("DimensionMismatch") # DimensionMismatch("dimensions must match: a has dims (Base.OneTo(2),), b has dims (Base.OneTo(3),), mismatch at 1")
    end

    c = ['I' for i in 1:length(a)]
    for i in 1:length(a)
        if a[i] == 'I'
            if b[i] == 'X'
                c[i] = 'X'
            elseif b[i] == 'Y'
                c[i] = 'Y'
            elseif b[i] == 'Z'
                c[i] = 'Z'
            end
        elseif a[i] == 'X'
            if b[i] == 'I'
                c[i] = 'X'
            elseif b[i] == 'Y'
                c[i] = 'Z'
            elseif b[i] == 'Z'
                c[i] = 'Y'
            end
        elseif a[i] == 'Y'
            if b[i] == 'I'
                c[i] = 'Y'
            elseif b[i] == 'X'
                c[i] = 'Z'
            elseif b[i] == 'Z'
                c[i] = 'X'
            end
        elseif a[i] == 'Z'
            if b[i] == 'I'
                c[i] = 'Z'
            elseif b[i] == 'X'
                c[i] = 'Y'
            elseif b[i] == 'Y'
                c[i] = 'X'
            end
        end
    end

    return c
end

# to get primitive element out of field in Nemo
function primitive_element(F::T; n_quo::Int = -1) where T <: Union{FqFiniteField, FqNmodFiniteField, Nemo.GaloisField, Nemo.GaloisFmpzField}
    n = order(F)-1
    k = fmpz(1)
    if n_quo != -1
        if !divisible(n, n_quo)
            return F(1)
        end
        n, k = ppio(n, fmpz(n_quo))
    end
    primefactors_n = collect(keys(factor(n).fac))

    x = rand(F)^Int(k)
    while iszero(x)
        x = rand(F)^Int(k)
    end
    while true
        found = true
        for l in primefactors_n
            if isone(x^Int(div(n, l)))
                found = false
                break
            end
        end
        if found
            break
        end
        x = rand(F)^Int(k)
        while iszero(x)
            x = rand(F)^Int(k)
        end
    end
    return x
end
