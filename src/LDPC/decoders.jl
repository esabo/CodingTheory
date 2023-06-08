# Example of using Gallager A and B (out should end up [1 1 1 0 0 0 0]
# H = matrix(GF(2), [1 1 0 1 1 0 0; 1 0 1 1 0 1 0; 0 1 1 1 0 0 1]);
# v = matrix(GF(2), 7, 1, [1, 1, 0, 0, 0, 0, 0]);
# out, iter, vtoc, ctov = messagepassing(H, v, missing, x->x, _GallagerAchecknodemessage, 100, :A);
# out, iter, vtoc, ctov = messagepassing(H, v, missing, x->x, _GallagerAchecknodemessage, 100, :B);

# TODO: scheduling
function messagepassing(H::T, v::T, chn::Union{Missing, MPNoiseModel},
    vtocmess::Function, ctovmess::Function,
    maxiter::Int=100, kind::Symbol=:SP, Bt::Int = 2) where T <: CTMatrixTypes

    kind ∈ (:SP, :MS, :A, :B) || throw(ArgumentError(":log or :prob expected for parameter kind"))
    kind ∈ (:SP, :MS) && chn == missing && throw(ArgumentError(":SP and :MS require a noise model"))
    Int(order(base_ring(H))) == Int(order(base_ring(v))) == 2 ||
        throw(ArgumentError("Currently only implemented for binary codes"))
    numcheck, numvar = size(H)
    numcheck > 0 && numvar > 0 || throw(ArgumentError("Input matrix of improper dimension"))
    (size(v) != (numvar, 1) && size(v) != (1, numvar)) && throw(ArgumentError("Vector has incorrect dimension"))
    
    HInt = FpmattoJulia(H)
    nrows(v) == 1 ? (w = FpmattoJulia(transpose(v));) : (w = FpmattoJulia(v);)
    checkadlist = [[] for _ in 1:numcheck]
    varadlist = [[] for _ in 1:numvar]
    maxiter += 1

    for r in 1:numcheck
        for c in 1:numvar
            if !iszero(HInt[r, c])
                push!(checkadlist[r], c)
                push!(varadlist[c], r)
            end
        end
    end

    S = kind ∈ (:A, :B) ? Int : Float64

    curr = zeros(Int, numvar)
    totals = zeros(S, 1, numvar)
    syn = zeros(Int, numcheck)
    checktovarmessages = zeros(S, numcheck, numvar, maxiter)
    vartocheckmessages = zeros(S, numvar, numcheck, maxiter)
    
    iter = 1
    if kind == :SP
        chninits = _channelinitBSC(w, chn.crossoverprob)
        if type == :BSC
            Threads.@threads for vn in 1:numvar
                vartocheckmessages[vn, varadlist[vn], 1] .= chninits[vn]
            end
        end
    elseif kind in (:A, :B)
        Threads.@threads for vn in 1:numvar
            vartocheckmessages[vn, varadlist[vn], :] .= w[vn]
        end
    end

    while iter < maxiter
        Threads.@threads for cn in 1:numcheck
            for v1 in checkadlist[cn]
                checktovarmessages[cn, v1, iter] = ctovmess(cn, v1, iter, checkadlist, vartocheckmessages)
            end
        end

        if kind == :SP
            Threads.@threads for vn in 1:numvar
                totals[vn] = chninit[vn]
                for c in varadlist[vn]
                    totals[vn] += checktovarmessages[c, vn, iter]
                end
            end
        end

        if kind == :SP
            # TODO: unroll or otherwise vectorize this operation
            @simd for i in 1:numvar
                curr[i] = total >= 0 ? 0 : 1
            end
        elseif kind in (:A, :B)
            @simd for i in 1:numvar
                len = length(varadlist[i])
                onecount = count(isone, view(checktovarmessages, varadlist[i], i, iter))
                d = fld(len, 2)
                curr[i] = onecount + (isone(w[i]) && iseven(len)) > d
            end
        end

        LinearAlgebra.mul!(syn, HInt, curr)
        @show curr
        @show syn .% 2
        iszero(syn .% 2) && return matrix(base_ring(H), 1, numvar, curr), iter, vartocheckmessages, checktovarmessages # others if necessary
        iter += 1

        if iter <= maxiter
            Threads.@threads for vn in 1:numvar
                for c1 in varadlist[vn]
                    # vartocheckmessages[vn, c1, iter] = vtocmess(vn, c1, chninit[vn], varadlist, checktovarmessages)
                    if kind == :SP
                        vartocheckmessages[vn, c1, iter] = totals[vn] - checktovarmessages[c1, vn, iter - 1]
                    elseif kind == :A && length(varadlist[vn]) > 1
                        if all(!isequal(w[vn]), checktovarmessages[c2, vn, iter - 1] for c2 in varadlist[vn] if c1 != c2)
                            vartocheckmessages[vn, c1, iter] ⊻= 1 
                        end
                    elseif kind == :B && length(varadlist[vn]) > Bt
                        if count(!isequal(w[vn]), checktovarmessages[c2, vn, iter - 1] for c2 in varadlist[vn] if c1 != c2) >= Bt
                            vartocheckmessages[vn, c1, iter] ⊻= 1 
                        end
                    end
                end
            end
        end
    end

    return # something
end

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

function _channelinitBSC(v::Vector{Int}, p::Float64)
    temp = log((1 - p) / p)
    chninit = zeros(Float64, 1, axes(v, 2))
    for i in axes(v, 2)
        chninit[i] = (-1)^v[i] * temp
    end
    return chninit
end

# made irrelevant by switching to the totals above
# function _SPvariablenodemessage(vn::Int, c1::Int, chninit::Float64, varadlist, checktovarmessages)
#     temp = chninit
#     for c2 in varadlist[vn]
#         if c2 != c1
#             temp += checktovarmessages[c2, vn, iter - 1]
#         end
#     end
#     return temp
# end

function _SPchecknodemessage(cn::Int, v1::Int, iter, checkadlist, vartocheckmessages)
    phi(x) = -log(tanh(0.5 * x))
    temp = 0.0
    s = 1
    for v2 in checkadlist[cn]
        if v2 != v1
            temp += phi(abs(vartocheckmessages[v2, cn, iter - 1]))
            s *= sign(vartocheckmessages[v2, cn, iter - 1])
        end
    end
    return s * phi(temp)
end

function _GallagerAchecknodemessage(cn::Int, v1::Int, iter::Int, checkadlist, vartocheckmessages)
    reduce(⊻, vartocheckmessages[v, cn, iter] for v in checkadlist[cn] if v != v1)
end
_GallagerBchecknodemessage(cn::Int, v1::Int, iter::Int, checkadlist, vartocheckmessages) = _GallagerAchecknodemessage(cn, v1, iter, checkadlist, vartocheckmessages)

# TODO: write functions for each kind of message
# BP
# Sum-product
# min-sum

# Mansour, Shanbhag, "Turbo Decoder Architectures for Low-Density Parity-Check Codes" (2002)
function messagepassingv2(H::T, v::T,
    vtocmess::Function, ctovmess::Function,
    maxiter::Int=100, kind::Symbol=:log) where T <: CTMatrixTypes

    kind ∈ (:log, :prob) || throw(ArgumentError(":log or :prob expected for parameter kind"))
    Int(order(base_ring(H))) == Int(order(base_ring(v))) == 2 ||
        throw(ArgumentError("Currently only implemented for binary codes"))
    numcheck, numvar = size(H)
    numcheck > 0 && numvar > 0 || throw(ArgumentError("Input matrix of improper dimension"))
    (size(v) != (numvar, 1) && size(v) != (1, numvar)) && throw(ArgumentError("Vector has incorrect dimension"))
    
    HInt = FpmattoJulia(H)
    nrows(v) == 1 ? (w = FpmattoJulia(transpose(v));) : (w = FpmattoJulia(v);)
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

    curr = zeros(Int, 1, numvar)
    syn = zeros(Int, 1, numcheck)
    checktovarmessages = zeros(Float64, numcheck, numvar, maxiter + 1)
    vartocheckmessages = zeros(Float64, numvar, maxiter + 1)
    temp = zeros(Float64, numvar, maxiter + 1)

    # initialize to channel
    @simd for i in 1:numvar
        temp[i, 1] = vartocheckmessages[i, 1] = ...
    end

    iter = 2
    while iter <= maxiter
        Threads.@threads for cn in 1:numcheck
            for v1 in checkadlist[cn]
                for v2 in checkadlist[cn]
                    if v2 != v1
                        checktovarmessages[cn, v1, iter] = ...
                    end
                end
                temp[v1, iter] = temp[v1, iter - 1] + checktovarmessages[cn, v1, iter]
            end
        end

        # if using probabilities, then > 1 implies evidence for hypothesis 1 over 2
        if kind == :log
            # TODO: unroll or otherwise vectorize this operation
            @simd for i in  1:numvar
                curr[i] = vartocheckmessages[i, ?, iter] >= 0 ? w[i] : (w[1] + 1) .% 2
            end
        else
            @simd for i in  1:numvar
                curr[i] = vartocheckmessages[i, ?, iter] >= 1 ? w[i] : (w[1] + 1) .% 2
            end
        end

        LinearAlgebra.mul!(syn, HInt, curr)
        iszero(syn .% 2) && return matrix(base_ring(H), 1, numvar, curr), iter # others if necessary

        vartocheckmessages = temp
        # initialize to channel
        @simd for i in 1:numvar
            temp[i] = ...
        end
        iter += 1
    end

    return # something
end

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
