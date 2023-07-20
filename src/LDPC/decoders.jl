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

function sumproduct(H::S, v::T, chn::MPNoiseModel, maxiter::Int=100) where {S <: CTMatrixTypes, T <: Union{Vector{<:Real}, CTMatrixTypes}}
    HInt, w, varadlist, checkadlist = _messagepassinginit(H, v, chn, maxiter, :SP, 2)
    return _messagepassing(HInt, w, chn, _SPchecknodemessage, varadlist, checkadlist, maxiter, :SP)
end

function sumproductboxplus(H::S, v::T, chn::MPNoiseModel, maxiter::Int=100) where {S <: CTMatrixTypes, T <: Union{Vector{<:Real}, CTMatrixTypes}}
    HInt, w, varadlist, checkadlist = _messagepassinginit(H, v, chn, maxiter, :SP, 2)
    return _messagepassing(HInt, w, chn, _SPchecknodemessageboxplus, varadlist, checkadlist, maxiter, :SP)
end

function minsum(H::S, v::T, chn::MPNoiseModel, maxiter::Int=100, attenuation::Float64 = 0.5) where {S <: CTMatrixTypes, T <: Vector{<:AbstractFloat}}
    HInt, w, varadlist, checkadlist = _messagepassinginit(H, v, chn, maxiter, :MS, 2)
    return _messagepassing(HInt, w, chn, _MSchecknodemessage, varadlist, checkadlist, maxiter, :MS, attenuation)
end

function _messagepassinginit(H::S, v::T, chn::Union{Missing, MPNoiseModel}, maxiter::Int, kind::Symbol, Bt::Int) where {S <: CTMatrixTypes, T <: Union{Vector{<:Real}, CTMatrixTypes}}
    kind ∈ (:SP, :MS, :A, :B) || throw(ArgumentError("Unknown value for parameter kind"))
    kind ∈ (:SP, :MS) && ismissing(chn) && throw(ArgumentError(":SP and :MS require a noise model"))
    Int(order(base_ring(H))) == 2 || throw(ArgumentError("Currently only implemented for binary codes"))
    numcheck, numvar = size(H)
    numcheck > 0 && numvar > 0 || throw(ArgumentError("Input matrix of improper dimension"))
    length(v) == numvar || throw(ArgumentError("Vector has incorrect dimension"))
    (kind == :B && !(1 <= Bt <= numcheck)) && throw(DomainError("Improper threshold for Gallager B"))
    2 <= maxiter || throw(DomainError("Number of maximum iterations must be at least two"))
    kind ∈ (:SP, :MS) && chn.type == :BAWGNC && !isa(v, Vector{<:AbstractFloat}) && throw(DomainError("Received message should be a vector of floats for BAWGNC."))
    kind ∈ (:SP, :MS) && chn.type == :BSC && !isa(v, Vector{Int}) && !isa(v, CTMatrixTypes) && throw(DomainError("Received message should be a vector of Ints for BSC."))
    
    HInt = FpmattoJulia(H)
    w = if T <: CTMatrixTypes
        Int.(data.(v)[:])
    else
        copy(v)
    end
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
function _messagepassing(H::Matrix{UInt64}, w::Vector{T}, chn::Union{Missing, MPNoiseModel},
    ctovmess::Function, varadlist::Vector{Vector{Any}}, checkadlist::Vector{Vector{Any}},
    maxiter::Int, kind::Symbol, Bt::Int=2, attenuation::Float64 = 0.5) where T <: Union{Int, AbstractFloat}

    numcheck, numvar = size(H)
    S = kind ∈ (:A, :B) ? Int : Float64
    curr = zeros(Int, numvar)
    if kind in (:SP, :MS)
        totals = zeros(S, 1, numvar)
    end
    syn = zeros(Int, numcheck)
    maxiter += 1 # probably should copy this
    checktovarmessages = zeros(S, numcheck, numvar, maxiter)
    vartocheckmessages = zeros(S, numvar, numcheck, maxiter)
    
    iter = 1
    if kind in (:SP, :MS)
        chninits = if chn.type == :BSC
            _channelinitBSC(w, chn.crossoverprob)
        elseif chn.type == :BAWGNC && kind == :SP
            _channelinitBAWGNCSP(w, chn.sigma)
        elseif chn.type == :BAWGNC && kind == :MS
            _channelinitBAWGNCMS(w)
        end
        for vn in 1:numvar
            vartocheckmessages[vn, varadlist[vn], 1] .= chninits[vn]
        end
    elseif kind in (:A, :B)
        for vn in 1:numvar
            vartocheckmessages[vn, varadlist[vn], :] .= w[vn]
        end
    end

    while iter < maxiter
        for cn in 1:numcheck
            for v1 in checkadlist[cn]
                checktovarmessages[cn, v1, iter] = ctovmess(cn, v1, iter, checkadlist, vartocheckmessages, attenuation)
            end
        end

        if kind in (:SP, :MS)
            for vn in 1:numvar
                totals[vn] = chninits[vn]
                for c in varadlist[vn]
                    totals[vn] += checktovarmessages[c, vn, iter]
                end
            end
        end

        if kind in (:SP, :MS)
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
                    if kind in (:SP, :MS)
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
    chninit = zeros(Float64, length(v))
    for i in eachindex(v)
        chninit[i] = (-1)^v[i] * temp
    end
    return chninit
end

function _channelinitBAWGNCSP(v::Vector{T}, σ::Float64) where T <: AbstractFloat
    temp = 2 / σ^2
    chninit = zeros(Float64, length(v))
    for i in eachindex(v)
        chninit[i] = temp * v[i]
    end
    return chninit
end

_channelinitBAWGNCMS(v::Vector{T}) where T <: AbstractFloat = v

function _SPchecknodemessage(cn::Int, v1::Int, iter, checkadlist, vartocheckmessages, atten = missing)
    phi(x) = -log(tanh(0.5 * x))
    temp = 0.0
    s = 1
    for v2 in checkadlist[cn]
        if v2 != v1
            x = vartocheckmessages[v2, cn, iter]
            # Note that x should never be 0 unless there is an erasure.
            # For now, this is coded as if there will never be an erasure.
            # This will simply error if x == 0.
            if x > 0
                temp += phi(x)
            else
                temp += phi(-x)
                s *= -1
            end
        end
    end
    return s * phi(temp)
end

⊞(a, b) = log((1 + exp(a + b)) / (exp(a) + exp(b)))
⊞(a...) = reduce(⊞, a...)
function _SPchecknodemessageboxplus(cn::Int, v1::Int, iter, checkadlist, vartocheckmessages, atten = missing)
    ⊞(vartocheckmessages[v2, cn, iter] for v2 in checkadlist[cn] if v2 != v1)
end

function _MSchecknodemessage(cn::Int, v1::Int, iter, checkadlist, vartocheckmessages, attenuation::Float64=0.5)
    temp = vartocheckmessages[checkadlist[cn][1], cn, iter]
    s = 1
    for v2 in checkadlist[cn]
        if v2 != v1
            x = vartocheckmessages[v2, cn, iter]
            # Note that x should never be 0
            if x > 0
                temp > x && (temp = x;)
            else
                temp > -x && (temp = -x;)
                s *= -1
            end
        end
    end
    return s * attenuation * temp
end

function _GallagerAchecknodemessage(cn::Int, v1::Int, iter::Int, checkadlist, vartocheckmessages, atten = missing)
    reduce(⊻, vartocheckmessages[v, cn, iter] for v in checkadlist[cn] if v != v1)
end
_GallagerBchecknodemessage(cn::Int, v1::Int, iter::Int, checkadlist, vartocheckmessages, atten = missing) = _GallagerAchecknodemessage(cn, v1, iter, checkadlist, vartocheckmessages, atten)

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

function _channeltoSNR(chn::MPNoiseModel)
    if cnh.type == :BAWGNC
        -10 * log10(chn.sigma^2)
    else
        throw(ArgumentError("Only supports BAWGNC currently."))
    end
end

function _channeltoSNR(type::Symbol, sigma::Real)
    if type == :BAWGNC
        -10 * log10(sigma^2)
    else
        throw(ArgumentError("Only supports BAWGNC currently."))
    end
end

function _initLPdecoderLDPC(H::Union{CTMatrixTypes, AbstractMatrix{<:Number}})
    checkadlist, _ = CodingTheory._nodeadjacencies(H)
    subsets = Vector{Vector{Vector{Int}}}()
    nr, nc = size(H)
    hasi = [[Vector{Int}() for _ in 1:nc] for _ in 1:nr]
    wmap = zeros(Int, nr)
    curr = 1
    for (j, cn) in enumerate(checkadlist)
        wmap[j] = curr
        innersubsets = Vector{Vector{Int}}()
        for S in powerset(cn)
            if iseven(length(S)) # && !isempty(S)
                push!(innersubsets, S)
                for i in S
                    push!(hasi[j][i], curr)
                end
                curr += 1
            end
        end
        push!(subsets, innersubsets)
    end

    model = Model(GLPK.Optimizer)
    @variable(model, 0 <= f[1:nc] <= 1)
    @variable(model, 0 <= w[1:curr - 1] <= 1)
    for i in 1:nr
        if i != nr
            @constraint(model, sum(w[wmap[i]:wmap[i + 1] - 1]) == 1)
        else
            @constraint(model, sum(w[wmap[i]:end]) == 1)
        end
    end
    for j in 1:nr
        for i in 1:nc
            if !isempty(hasi[j][i])
                @constraint(model, f[i] == sum(w[hasi[j][i]]))
            end
        end
    end
    @objective(model, Min, sum(0 * f[i] for i in 1:nc))
    return model
end
_initLPdecoderLDPC(C::AbstractLinearCode) = _initLPdecoderLDPC(paritycheckmatrix(C))

function _LPdecoderLDPC(model::JuMP.Model, v::Union{CTMatrixTypes, Vector{<:Integer}}, Ch::CodingTheory.BinarySymmetricChannel)
    γ = CodingTheory._channelinitBSC(isa(v, Vector) ? v : Int.(data.(v))[:], Ch.param)
    @objective(model, Min, dot(γ, model[:f]))
    optimize!(model)
    termination_status(model) == MOI.INFEASIBLE && throw(DomainError("No solution exists"))
    @assert termination_status(model) == MOI.OPTIMAL "Didn't find an optimal point"
    w = value.(model[:f])
    all(isinteger(x) for x in w) || @warn "Solution is not integral"
    return w
end

function LPdecoderLDPC(H::Union{CTMatrixTypes, AbstractMatrix{<:Number}},
                       v::Union{CTMatrixTypes, Vector{<:Integer}},
                       Ch::BinarySymmetricChannel)
    model = _initLPdecoderLDPC(H)
    return _LPdecoderLDPC(model, v, Ch)
end
LPdecoderLDPC(C::AbstractLinearCode, v::Union{CTMatrixTypes, Vector{<:Integer}}, Ch::CodingTheory.BinarySymmetricChannel) = LPdecoderLDPC(paritycheckmatrix(C), v, Ch)


function decodersimulation(H::CTMatrixTypes, decoder::Symbol, noisetype::Symbol,
                           noise::Union{Vector{T}, AbstractRange{T}} where T<:Real,
                           maxiter::Int=100, numruns::Int=100000, seed::Int = 123,
                           attenuation::Float64 = 0.5)

    decodersimulation(H, decoder, noisetype, noise, maxiter,
                      [numruns for n in noise], seed, attenuation)
end

function decodersimulation(H::CTMatrixTypes, decoder::Symbol, noisetype::Symbol,
                           noise::Union{Vector{T}, AbstractRange{T}} where T<:Real,
                           maxiter::Int=100, numruns::Vector{Int} = [100000 for n in noise],
                           seed::Int = 123, attenuation::Float64 = 0.5)

    decoder in (:A, :B, :SP, :MS) || throw(ArgumentError("Unsupported decoder"))
    noisetype in (:BSC, :BAWGNC) || throw(ArgumentError("Only supports BSC and BAWGNC"))
    decoder in (:A, :B) && noisetype == :BAWGNC && throw(ArgumentError("BAWGNC not supported for Gallager decoders."))
    0 <= minimum(noise) || throw(ArgumentError("Must have non-negative noise"))
    maximum(noise) > 1 && noisetype == :BSC && throw(ArgumentError("Crossover probability must be in the range [0,1]"))

    # We use an explicit pseudoRNG with the given seed.
    # Note that Xoshiro is default in Julia.
    # Also note that threading breaks the reproducibility.
    rng = Xoshiro(seed)

    FER = zeros(length(noise))
    BER = zeros(length(noise))
    n = ncols(H)
    cnmsg = decoder == :SP ? _SPchecknodemessage : _MSchecknodemessage

    for k in eachindex(noise)
        chn = MPNoiseModel(noisetype, noise[k])
        w = noisetype == :BSC ? zeros(Int, n) : ones(n)
        HInt, _, varadlist, checkadlist = _messagepassinginit(H, w, chn, maxiter, decoder, 2)
        FEtotal = 0 # number of frame errors
        BEtotal = 0 # number of bit errors
        @threads for i in 1:numruns[k]
        # for i in 1:numruns
            for j in 1:n
                if noisetype == :BSC
                    w[j] = Int(rand(rng) < chn.crossoverprob)
                else # BAWGNC
                    w[j] = 1.0 + randn(rng, Float64) * chn.sigma
                end
            end
            iszero(w) && continue
            flag, curr, _, _, _ = _messagepassing(HInt, w, chn, cnmsg, varadlist, checkadlist, maxiter, decoder, 2, attenuation)
            if !(flag && iszero(curr))
                FEtotal += 1
                BEtotal += count(!iszero, curr)
            end
        end
        FER[k] = FEtotal / numruns[k]
        BER[k] = BEtotal / (numruns[k] * n)
    end
    return FER, BER
end

function decodersimulation2(H::CTMatrixTypes, decoder::Symbol, noisetype::Symbol,
                            noise::Union{Vector{T}, AbstractRange{T}} where T<:Real,
                            maxiter::Int=100, numruns::Int = 1000,
                            seed::Int = 123, attenuation::Float64 = 1.0)

    decoder in (:A, :B, :SP, :MS) || throw(ArgumentError("Unsupported decoder"))
    noisetype == :BSC || throw(ArgumentError("Only supports BSC"))
    0 <= minimum(noise) || throw(ArgumentError("Must have non-negative noise"))
    maximum(noise) > 1 && noisetype == :BSC && throw(ArgumentError("Crossover probability must be in the range [0,1]"))

    # we'll use an explicit pseudoRNG with the given seed. Note that Xoshiro is default in Julia.
    rng = Xoshiro(seed)

    FER = zeros(length(noise))
    BER = zeros(length(noise))
    ε = zeros(length(noise))
    n = ncols(H)
    cnmsg = decoder == :SP ? _SPchecknodemessage : _MSchecknodemessage

    for k in eachindex(noise)

        # p[i] is the probability of having i - 1 bit errors
        temp = BigFloat(noise[k])
        p = BigFloat[temp^i * (1 - temp)^(n - i) * binomial(big(n), big(i)) /
                     sum(temp^j * (1 - temp)^(n - j) * binomial(big(n), big(j))
                         for j in 0:n) for i in 0:n]
        p_partialsum = [sum(p[j] for j in 1:i) for i in 1:length(p)]
        maxnerr = max(findfirst(p_partialsum .>= 1 - BigFloat("1e-7")) - 1, 6)
        # @show maxnerr
        ε[k] = 1 - p_partialsum[maxnerr + 1]

        chn = MPNoiseModel(noisetype, noise[k])
        w = noisetype == :BSC ? zeros(Int, n) : ones(n)
        HInt, _, varadlist, checkadlist = _messagepassinginit(H, w, chn, maxiter, decoder, 2)

        FEtotal = zeros(Int, maxnerr)
        BEtotal = zeros(Int, maxnerr)
        FER[k] = p[1]
        BER[k] = p[1]

        for e in 1:maxnerr

            # importance sampling:
            numruns_for_e = Int(cld(numruns * p[e+1], sum(p[i] for i in 2:maxnerr+1)))

            # naive:
            # numruns_for_e = Int(cld(numruns, maxnerr))

            for i in 1:numruns_for_e
                w = zeros(Int, n)
                w[shuffle(rng, 1:n)[1:e]] .= 1
                flag, curr, _, _, _ = _messagepassing(HInt, w, chn, cnmsg, varadlist, checkadlist, maxiter, decoder, 2, attenuation)
                if !(flag && iszero(curr))
                    FEtotal[e] += 1
                    BEtotal[e] += count(!iszero, curr)
                end
            end
            FER[k] += p[e+1] * (numruns_for_e - FEtotal[e]) / numruns_for_e
            BER[k] += p[e+1] * (numruns_for_e * n - BEtotal[e]) / (numruns_for_e * n)
        end
        FER[k] = 1 - FER[k]
        BER[k] = 1 - BER[k]
    end
    return FER, BER, ε
end

using Plots: plot, savefig
function testsimulation()
    # C = BCHCode(2, 2^10-1, 53, 1);
    # C = BCHCode(2, 103, 3);
    # H = paritycheckmatrix(C);

    H = matrix(GF(2), 10, 20, [1 0 1 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0;
                               0 1 0 1 0 1 1 0 0 0 0 1 0 0 0 0 0 0 0 0;
                               0 0 1 0 1 0 1 1 0 0 0 0 1 0 0 0 0 0 0 0;
                               1 0 0 1 0 0 0 1 1 0 0 0 0 1 0 0 0 0 0 0;
                               0 1 0 0 1 0 0 0 1 1 0 0 0 0 1 0 0 0 0 0;
                               1 1 0 0 0 0 1 0 1 0 0 0 0 0 0 1 0 0 0 0;
                               0 1 1 0 0 0 0 1 0 1 0 0 0 0 0 0 1 0 0 0;
                               0 0 1 1 0 1 0 0 1 0 0 0 0 0 0 0 0 1 0 0;
                               0 0 0 1 1 0 1 0 0 1 0 0 0 0 0 0 0 0 1 0;
                               1 0 0 0 1 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1]);

    p0 = 0.0001:0.0001:0.0009;
    p1 = 0.001:0.001:0.01;
    p2 = 0.025:0.025:0.15;
    @time FER0, BER0 = decodersimulation(H, :SP, :BSC, p0, 100, 20000000, 123);
    @time FER1, BER1 = decodersimulation(H, :SP, :BSC, p1, 100, 500000, 123);
    @time FER2, BER2 = decodersimulation(H, :SP, :BSC, p2, 100, 10000, 123);
    @time FER0ms, BER0ms = decodersimulation(H, :MS, :BSC, p0, 100, 20000000, 123);
    @time FER1ms, BER1ms = decodersimulation(H, :MS, :BSC, p1, 100, 500000, 123);
    @time FER2ms, BER2ms = decodersimulation(H, :MS, :BSC, p2, 100, 10000, 123);
    @time FER0msnoa, BER0msnoa = decodersimulation(H, :MS, :BSC, p0, 100, 20000000, 123, 1.);
    @time FER1msnoa, BER1msnoa = decodersimulation(H, :MS, :BSC, p1, 100, 500000, 123, 1.);
    @time FER2msnoa, BER2msnoa = decodersimulation(H, :MS, :BSC, p2, 100, 10000, 123, 1.);
    p = vcat(p0, p1, p2);
    FER = vcat(FER0, FER1, FER2);
    BER = vcat(BER0, BER1, BER2);
    FERms = vcat(FER0ms, FER1ms, FER2ms);
    BERms = vcat(BER0ms, BER1ms, BER2ms);
    FERmsnoa = vcat(FER0msnoa, FER1msnoa, FER2msnoa);
    BERmsnoa = vcat(BER0msnoa, BER1msnoa, BER2msnoa);

    p = 10 .^ collect(-4:.2:-0.8);
    @time F1, _ = CodingTheory.decodersimulation2(H, :SP, :BSC, p, 100, 10000, 1, 1.0);
    @time F2, _ = CodingTheory.decodersimulation2(H, :MS, :BSC, p, 100, 5000, 1, 1.0);
    @time F3, _ = CodingTheory.decodersimulation2(H, :MS, :BSC, p, 100, 5000, 1, 0.5);
    @time F4, _ = CodingTheory.decodersimulation2(H, :MS, :BSC, p, 100, 5000, 1, 0.1);

    plt = plot(log10.(p), log10.([F1 F2 F3 F4]),
               label = ["FER, SP" "FER, MS atten=1.0" "FER, MS atten=0.5" "FER, MS atten=0.1"],
               xlabel = "Crossover probability",
               ylabel = "Error rate",
               title = "BSC with a [20,10,5] code",
               # xlims = (0, maximum(p) * 1.02),
               # ylims = (0, max(maximum(FER), maximum(BER)) * 1.02),
               # xticks = (-4:-1, ["1e-4", "1e-3", "1e-2", "1e-1"]),
               # yticks = (-6:0, ["1e-6", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1", "1e0"]),
               # yscale = :log,
               marker = :dot);
    xticks!(plt, (xticks(plt)[1][1], "1e" .* xticks(plt)[1][2]));
    yticks!(plt, (yticks(plt)[1][1], "1e" .* yticks(plt)[1][2]));


    # σ = 0.1:0.1:1
    # FER, BER = decodersimulation(H, :SP, :BAWGNC, p, 100, 100, 123);
    # SNR = CodingTheory._channeltoSNR.(:BAWGNC, σ)
    # plt = plot(SNR, [FER BER],
    #            label = ["FER" "BER"],
    #            xlabel = "Noise (dB)",
    #            ylabel = "Error rate",
    #            marker = :dot);


    savefig(plt, "test.png");
end


# function profiletest()
#     C = BCHCode(2,103,3);
#     H = paritycheckmatrix(C);
#     chn = MPNoiseModel(:BSC, 0.01);
#     v = zero_matrix(C.F, 1, 103);
#     v[1,1] = 1;
#     sumproduct(H, v, chn);
#     Profile.clear()
#     @profile sumproduct(H, v, chn)
# end
