# Example of using Gallager A and B (out should end up [1 1 1 0 0 0 0]
# H = matrix(GF(2), [1 1 0 1 1 0 0; 1 0 1 1 0 1 0; 0 1 1 1 0 0 1]);
# v = matrix(GF(2), 7, 1, [1, 1, 0, 0, 0, 0, 0]);
# flag, out, iter, vtoc, ctov = GallagerA(H, v, 100);
# flag, out, iter, vtoc, ctov = GallagerB(H, v, 100);
# nm = MPNoiseModel(:BSC, 1/7)
# flag, out, iter, vtoc, ctov = sumproduct(H, v, nm, 100);

struct MPNoiseModel
    type::Symbol
    crossoverprob::Union{Float64, Missing}
    sigma::Union{Float64, Missing}
end

function MPNoiseModel(type::Symbol, x::Float64)
    if type == :BSC
        return MPNoiseModel(type, x, missing)
    elseif type == :BAWGNC
        return MPNoiseModel(type, missing, x)
    else
        throw(ArgumentError("Unsupported noise model type $type"))
    end
end

function GallagerA(H::T, v::T, maxiter::Int=100) where T <: CTMatrixTypes
    HInt, w, varadlist, checkadlist = _messagepassinginit(H, v, missing, maxiter, :A, 2)
    return _messagepassing(HInt, w, missing, _GallagerAchecknodemessage, varadlist, checkadlist, maxiter, :A)
end

function GallagerB(H::T, v::T, maxiter::Int=100, threshold::Int=2) where T <: CTMatrixTypes
    HInt, w, varadlist, checkadlist = _messagepassinginit(H, v, missing, maxiter, :B, threshold)
    return _messagepassing(HInt, w, missing, _GallagerBchecknodemessage, varadlist, checkadlist, maxiter, :B, threshold)
end

function sumproduct(H::T, v::T, chn::MPNoiseModel, maxiter::Int=100) where T <: CTMatrixTypes
    HInt, w, varadlist, checkadlist = _messagepassinginit(H, v, chn, maxiter, :SP, 2)
    return _messagepassing(HInt, w, chn, _SPchecknodemessage, varadlist, checkadlist, maxiter, :SP)
end

function _messagepassinginit(H::T, v::T, chn::Union{Missing, MPNoiseModel}, maxiter::Int, kind::Symbol, Bt::Int) where T <: CTMatrixTypes
    kind ∈ (:SP, :MS, :A, :B) || throw(ArgumentError("Unknown value for parameter kind"))
    kind ∈ (:SP, :MS) && ismissing(chn) && throw(ArgumentError(":SP and :MS require a noise model"))
    Int(order(base_ring(H))) == Int(order(base_ring(v))) == 2 ||
        throw(ArgumentError("Currently only implemented for binary codes"))
    numcheck, numvar = size(H)
    numcheck > 0 && numvar > 0 || throw(ArgumentError("Input matrix of improper dimension"))
    (size(v) != (numvar, 1) && size(v) != (1, numvar)) && throw(ArgumentError("Vector has incorrect dimension"))
    (kind == :B && !(1 <= Bt <= numcheck)) && throw(DomainError("Improper threshold for Gallager B"))
    2 <= maxiter || throw(DomainError("Number of maximum iterations must be at least two"))
    
    HInt = FpmattoJulia(H)
    w = Int.(data.(v)[:])
    checkadlist = [[] for _ in 1:numcheck]
    varadlist = [[] for _ in 1:numvar]

    for r in 1:numcheck
        for c in 1:numvar
            if !iszero(HInt[r, c])
                push!(checkadlist[r], c)
                push!(varadlist[c], r)
            end
        end
    end

    return HInt, w, varadlist, checkadlist
end

# TODO: scheduling
function _messagepassing(H::Matrix{UInt64}, w::Vector{Int}, chn::Union{Missing, MPNoiseModel}, ctovmess::Function,
    varadlist::Vector{Vector{Any}}, checkadlist::Vector{Vector{Any}}, maxiter::Int, kind::Symbol, Bt::Int=2)

    numcheck, numvar = size(H)
    S = kind ∈ (:A, :B) ? Int : Float64
    curr = zeros(Int, numvar)
    if kind == :SP
        totals = zeros(S, 1, numvar)
    end
    syn = zeros(Int, numcheck)
    maxiter += 1 # probably should copy this
    checktovarmessages = zeros(S, numcheck, numvar, maxiter)
    vartocheckmessages = zeros(S, numvar, numcheck, maxiter)
    
    iter = 1
    if kind == :SP
        if chn.type == :BSC
            chninits = _channelinitBSC(w, chn.crossoverprob)
            for vn in 1:numvar
                vartocheckmessages[vn, varadlist[vn], 1] .= chninits[vn]
            end
        end
    elseif kind in (:A, :B)
        for vn in 1:numvar
            vartocheckmessages[vn, varadlist[vn], :] .= w[vn]
        end
    end

    while iter < maxiter
        for cn in 1:numcheck
            for v1 in checkadlist[cn]
                checktovarmessages[cn, v1, iter] = ctovmess(cn, v1, iter, checkadlist, vartocheckmessages)
            end
        end

        if kind == :SP
            for vn in 1:numvar
                totals[vn] = chninits[vn]
                for c in varadlist[vn]
                    totals[vn] += checktovarmessages[c, vn, iter]
                end
            end
        end

        if kind == :SP
            @simd for i in 1:numvar
                curr[i] = totals[i] >= 0 ? 0 : 1
            end
        elseif kind in (:A, :B)
            @simd for i in 1:numvar
                len = length(varadlist[i])
                onecount = count(isone, view(checktovarmessages, varadlist[i], i, iter))
                d = fld(len, 2)
                curr[i] = onecount + (isone(w[i]) && iseven(len)) > d
            end
        end

        LinearAlgebra.mul!(syn, H, curr)
        # @show curr
        # @show syn .% 2
        iszero(syn .% 2) && return true, curr, iter, vartocheckmessages, checktovarmessages # others if necessary
        iter += 1

        if iter <= maxiter
            for vn in 1:numvar
                for c1 in varadlist[vn]
                    if kind == :SP
                        vartocheckmessages[vn, c1, iter] = totals[vn] - checktovarmessages[c1, vn, iter - 1]
                    elseif kind == :A && length(varadlist[vn]) > 1
                        if all(!Base.isequal(w[vn]), checktovarmessages[c2, vn, iter - 1] for c2 in varadlist[vn] if c1 != c2)
                            vartocheckmessages[vn, c1, iter] ⊻= 1 
                        end
                    elseif kind == :B && length(varadlist[vn]) >= Bt
                        if count(!Base.isequal(w[vn]), checktovarmessages[c2, vn, iter - 1] for c2 in varadlist[vn] if c1 != c2) >= Bt
                            vartocheckmessages[vn, c1, iter] ⊻= 1 
                        end
                    end
                end
            end
        end
    end

    return false, curr, iter, vartocheckmessages, checktovarmessages
end

function _channelinitBSC(v::Vector{T}, p::Float64) where T <: Integer
    temp = log((1 - p) / p)
    chninit = zeros(Float64, length(v), 1)
    for i in eachindex(v)
        chninit[i] = (-1)^v[i] * temp
    end
    return chninit
end

function _SPchecknodemessage(cn::Int, v1::Int, iter, checkadlist, vartocheckmessages)
    phi(x) = -log(tanh(0.5 * x))
    temp = 0.0
    s = 1
    for v2 in checkadlist[cn]
        if v2 != v1
            temp += phi(abs(vartocheckmessages[v2, cn, iter]))
            s *= sign(vartocheckmessages[v2, cn, iter])
        end
    end
    return s * phi(temp)
end

function _GallagerAchecknodemessage(cn::Int, v1::Int, iter::Int, checkadlist, vartocheckmessages)
    reduce(⊻, vartocheckmessages[v, cn, iter] for v in checkadlist[cn] if v != v1)
end
_GallagerBchecknodemessage(cn::Int, v1::Int, iter::Int, checkadlist, vartocheckmessages) = _GallagerAchecknodemessage(cn, v1, iter, checkadlist, vartocheckmessages)

# min-sum

# Mansour, Shanbhag, "Turbo Decoder Architectures for Low-Density Parity-Check Codes" (2002)


function findMPschedule(H::CodingTheory.CTMatrixTypes)
    numcheck, numvar = size(H)
    numcheck > 0 && numvar > 0 || throw(ArgumentError("Input matrix of improper dimension"))

    checkadlist = [[] for _ in 1:numcheck]
    for r in 1:numcheck
        for c in 1:numvar
            iszero(H[r, c]) || push!(checkadlist[r], c)
        end
    end

    schedlist = [[1]]
    for cn in 2:numcheck
        found = false
        for sched in schedlist
            if !any(x ∈ checkadlist[y] for y in sched for x ∈ checkadlist[cn])
                push!(sched, cn)
                sort!(schedlist, lt=(x, y) -> length(x) < length(y))
                found = true
                break
            end
        end
        !found && push!(schedlist, [cn])
    end
    
    return schedlist
end

# # need to take in a range of parameter values here
function sumproductsimulation(H::CTMatrixTypes, chn::MPNoiseModel, maxiter::Int=100, numruns::Int=100000)
    # check on parameters here
    HInt, _, varadlist, checkadlist = _messagepassinginit(H, zero_matrix(base_ring(H), 1, ncols(H)), chn, maxiter, :SP, 2)
    failure = 0
    # sample here for numruns copies of w
    n = ncols(H)
    w = zeros(Int, n)
    for i in 1:numruns
        for j in 1:n
            w[j] = Int(rand() < chn.crossoverprob)
        end
        iszero(w) && continue
        flag, curr, _, _, _ = _messagepassing(HInt, w, chn, _SPchecknodemessage, varadlist, checkadlist, maxiter, :SP)
        (flag && iszero(curr)) || (failure += 1;)
    end
    return failure / numruns
end

# function simulation(type::Symbol = :SP, H::CTMatrixTypes, parr::Vector{Float64}, maxiter::Int = 100, numrus::Int = 100000)
#     if type == :SP
#         for p in parr
#             chn = MPNoiseModel(:BSC, p)
#         end
#     elseif type == :A
#     end
# end
