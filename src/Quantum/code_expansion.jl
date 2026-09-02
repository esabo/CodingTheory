

"""
    laplacian_matrix(G::Grphs.SimpleGraph)

Return the Laplacian matrix of the graph `G` (Degree - Adjacency).
"""
function laplacian_matrix(G::Grphs.SimpleGraph)
    A = Grphs.adjacency_matrix(G)
    D = spdiagm(0 => Grphs.degree(G))
    return D - A
end

"""
    spectral_gap(G::Grphs.SimpleGraph)

Return the second-largest eigenvalue (in absolute value) of the graph `G`.
"""
function spectral_gap(G::Grphs.SimpleGraph)
    A = Grphs.adjacency_matrix(G)
    evals = eigvals(Symmetric(Matrix(A)))
    
    max_eval = maximum(abs.(evals))
    non_trivial_evals = filter(x -> abs(x) < max_eval - 1e-7, evals)
    
    isempty(non_trivial_evals) && return 0.0
    return maximum(abs.(non_trivial_evals))
end

"""
    spectral_gap(S::AbstractSubsystemCode, check_type::Symbol=:both)

Return the spectral gap of the code's Tanner graph.
For CSS codes, `check_type` can be `:X`, `:Z`, or `:both` to analyze the respective sub-graphs.
Results are cached.
"""
function spectral_gap(S::AbstractSubsystemCode, check_type::Symbol=:both)
    if CSSTrait(typeof(S)) == IsNotCSS()
        check_type == :both || @warn "check_type ignored for non-CSS codes. Using full graph."
        haskey(S.cache, :spectral_gap) && return S.cache[:spectral_gap]
        
        # Extract the graph from the cached tuple
        G = Tanner_graph(S)[1]
        gap = spectral_gap(G)
        S.cache[:spectral_gap] = gap
        return gap
    end

    check_type ∈ (:X, :Z, :both) || throw(ArgumentError("check_type must be :X, :Z, or :both"))
    
    cache_key = Symbol("spectral_gap_", check_type)
    haskey(S.cache, cache_key) && return S.cache[cache_key]
    
    if check_type == :X
        G = Tanner_graph_X(S)[1]
    elseif check_type == :Z
        G = Tanner_graph_Z(S)[1]
    else
        G = Tanner_graph(S)[1]
    end
    
    gap = spectral_gap(G)
    S.cache[cache_key] = gap
    return gap
end

"""
    girth(S::AbstractSubsystemCode, check_type::Symbol=:both)

Return the girth (shortest cycle) of the code's Tanner graph.
For CSS codes, `check_type` can be `:X`, `:Z`, or `:both` to analyze the respective sub-graphs.
Results are cached.
"""
function girth(S::AbstractSubsystemCode, check_type::Symbol=:both)
    if CSSTrait(typeof(S)) == IsNotCSS()
        check_type == :both || @warn "check_type ignored for non-CSS codes. Using full graph."
        haskey(S.cache, :girth) && return S.cache[:girth]
        
        G = Tanner_graph(S)[1]
        g = Grphs.girth(G)
        S.cache[:girth] = g
        return g
    end
    
    check_type ∈ (:X, :Z, :both) || throw(ArgumentError("check_type must be :X, :Z, or :both"))
    
    cache_key = Symbol("girth_", check_type)
    haskey(S.cache, cache_key) && return S.cache[cache_key]
    
    if check_type == :X
        G = Tanner_graph_X(S)[1]
    elseif check_type == :Z
        G = Tanner_graph_Z(S)[1]
    else
        G = Tanner_graph(S)[1]
    end
    
    g = Grphs.girth(G)
    S.cache[cache_key] = g
    return g
end

"""
    estimated_QLTC_soundness(S::AbstractCSSCode, check_type::Symbol=:X; max_wt::Int=3, samples::Int=1000)

Estimates the quantum soundness (coboundary expansion) of the code by sampling random errors 
up to `max_wt` and calculating the minimum ratio of syndrome weight to error weight.
Since exact soundness is NP-hard, this provides a heuristic upper bound.
"""
function estimated_QLTC_soundness(S::AbstractCSSCode, check_type::Symbol=:X; max_wt::Int=3, samples::Int=1000)
    check_type ∈ (:X, :Z) || throw(ArgumentError("check_type must be :X or :Z"))
    
    H = check_type == :X ? X_stabilizers(S) : Z_stabilizers(S)
    n = length(S)
    F = field(S)
    
    min_soundness = Inf
    
    for wt in 1:max_wt
        for _ in 1:samples
            # Generate a random error of specific weight
            E = zero_matrix(F, n, 1)
            indices = randperm(n)[1:wt]
            for idx in indices
                E[idx, 1] = F(1)
            end
            
            # Calculate syndrome weight
            syndrome = H * E
            syn_wt = count(!iszero, syndrome)
            
            ratio = syn_wt / wt
            if ratio < min_soundness
                min_soundness = ratio
            end
        end
    end
    
    return min_soundness
end

"""
    spectral_soundness_bound(S::AbstractCSSCode, check_type::Symbol=:X)

Returns a theoretical lower bound on the soundness of the code based on the spectral gap 
of its Tanner graph. 
"""
function spectral_soundness_bound(S::AbstractCSSCode, check_type::Symbol=:X)
    gap = spectral_gap(S, check_type)
    
    # Extract max check degree
    H = check_type == :X ? X_stabilizers(S) : Z_stabilizers(S)
    d_c = maximum(count(!iszero, H[i, :]) for i in 1:size(H, 1))
    
    # A standard spectral bound for expansion (Alon-Boppana / Cheeger variants)
    # The exact formula depends on the regularity of the graph.
    # For a d-regular graph, expansion is bounded by (d - λ) / 2
    return (d_c - gap) / 2
end

using Base.Threads: @threads, nthreads, threadid, Atomic, atomic_xchg!

# -----------------------------------------------------------------------------
# Bit-Packing Helper
# -----------------------------------------------------------------------------
function _bitpack_matrix(H::CTMatrixTypes)
    nr, nc = size(H)
    num_words = cld(nr, 64)
    H_packed = zeros(UInt64, nc, num_words)
    
    for c in 1:nc
        for r in 1:nr
            if !iszero(H[r, c])
                word_idx = (r - 1) ÷ 64 + 1
                bit_idx = (r - 1) % 64
                H_packed[c, word_idx] |= (UInt64(1) << bit_idx)
            end
        end
    end
    return H_packed, num_words, nr, nc
end

# -----------------------------------------------------------------------------
# 1. Deterministic Soundness
# -----------------------------------------------------------------------------
function _QLTC_soundness_dfs!(start_col::Int, depth::Int, max_wt::Int, nc::Int, num_words::Int, 
                         Δ::Int, syndromes::Matrix{UInt64}, H_packed::Matrix{UInt64}, 
                         min_ratios::Vector{Float64}, tid::Int)
    if depth == max_wt
        return
    end

    syn_wt = 0
    @inbounds for w in 1:num_words
        syn_wt += count_ones(syndromes[depth + 1, w])
    end

    best_future_ratio = Inf
    for k in (depth + 1):max_wt
        best_possible_syn_wt = max(0, syn_wt - (k - depth) * Δ)
        ratio_bound = best_possible_syn_wt / k
        if ratio_bound < best_future_ratio
            best_future_ratio = ratio_bound
        end
    end

    @inbounds if best_future_ratio >= min_ratios[tid]
        return
    end

    for c in start_col:nc
        next_syn_wt = 0
        @inbounds for w in 1:num_words
            val = syndromes[depth + 1, w] ⊻ H_packed[c, w]
            syndromes[depth + 2, w] = val
            next_syn_wt += count_ones(val)
        end

        if next_syn_wt > 0
            ratio = next_syn_wt / (depth + 1)
            @inbounds if ratio < min_ratios[tid]
                min_ratios[tid] = ratio
            end
        end

        _QLTC_soundness_dfs!(c + 1, depth + 1, max_wt, nc, num_words, Δ, syndromes, H_packed, min_ratios, tid)
    end
end

"""
    deterministic_QLTC_soundness(S::AbstractSubsystemCode, check_type::Symbol=:X; max_wt::Int=4, upper_bound::Float64=Inf)

Returns the exact worst-case soundness ratio (syndrome weight / error weight) for 
all errors up to `max_wt`. You can supply an `upper_bound` to instantly trigger deep Branch and Bound pruning.
"""
function deterministic_QLTC_soundness(S::AbstractSubsystemCode, check_type::Symbol=:X; max_wt::Int=4, upper_bound::Float64=Inf)
    check_type ∈ (:X, :Z) || throw(ArgumentError("check_type must be :X or :Z"))
    H_orig = check_type == :X ? X_stabilizers(S) : Z_stabilizers(S)
    
    col_wts = [count(!iszero, H_orig[:, c]) for c in 1:size(H_orig, 2)]
    p = sortperm(col_wts)
    H = H_orig[:, p]
    Δ = maximum(col_wts)
    
    H_packed, num_words, nr, nc = _bitpack_matrix(H)
    
    n_threads = nthreads()
    min_ratios = fill(upper_bound, n_threads)
    thread_syndromes = [zeros(UInt64, max_wt + 1, num_words) for _ in 1:n_threads]

    @threads for c1 in 1:nc
        tid = threadid()
        syndromes = thread_syndromes[tid]
        
        syn_wt = 0
        @inbounds for w in 1:num_words
            val = H_packed[c1, w]
            syndromes[2, w] = val
            syn_wt += count_ones(val)
        end
        
        if syn_wt > 0
            ratio = syn_wt / 1.0
            @inbounds if ratio < min_ratios[tid]
                min_ratios[tid] = ratio
            end
        end
        
        _QLTC_soundness_dfs!(c1 + 1, 1, max_wt, nc, num_words, Δ, syndromes, H_packed, min_ratios, tid)
    end
    
    return minimum(min_ratios)
end

# -----------------------------------------------------------------------------
# 2. Confinement Profile
# -----------------------------------------------------------------------------
function _confinement_dfs!(start_col::Int, depth::Int, max_wt::Int, nc::Int, num_words::Int, 
                           Δ::Int, syndromes::Matrix{UInt64}, H_packed::Matrix{UInt64}, 
                           profile::Vector{Int}, max_syndrome_wt::Int)
    if depth == max_wt
        return
    end

    syn_wt = 0
    @inbounds for w in 1:num_words
        syn_wt += count_ones(syndromes[depth + 1, w])
    end

    best_possible_future_syn_wt = max(0, syn_wt - (max_wt - depth) * Δ)
    if best_possible_future_syn_wt > max_syndrome_wt
        return
    end

    for c in start_col:nc
        next_syn_wt = 0
        @inbounds for w in 1:num_words
            val = syndromes[depth + 1, w] ⊻ H_packed[c, w]
            syndromes[depth + 2, w] = val
            next_syn_wt += count_ones(val)
        end

        if next_syn_wt <= max_syndrome_wt
            @inbounds if depth + 1 < profile[next_syn_wt + 1]
                profile[next_syn_wt + 1] = depth + 1
            end
        end

        _confinement_dfs!(c + 1, depth + 1, max_wt, nc, num_words, Δ, syndromes, H_packed, profile, max_syndrome_wt)
    end
end

"""
    confinement_profile(S::AbstractSubsystemCode, check_type::Symbol=:X; max_wt::Int=4, max_syndrome_wt::Int=typemax(Int))

Returns a deterministic dictionary mapping each observed syndrome weight to the 
minimum error weight (up to `max_wt`) that can trigger it. Supplying `max_syndrome_wt` heavily prunes the search space.
"""
function confinement_profile(S::AbstractSubsystemCode, check_type::Symbol=:X; max_wt::Int=4, max_syndrome_wt::Int=typemax(Int))
    check_type ∈ (:X, :Z) || throw(ArgumentError("check_type must be :X or :Z"))
    H_orig = check_type == :X ? X_stabilizers(S) : Z_stabilizers(S)
    
    col_wts = [count(!iszero, H_orig[:, c]) for c in 1:size(H_orig, 2)]
    p = sortperm(col_wts)
    H = H_orig[:, p]
    Δ = maximum(col_wts)
    
    H_packed, num_words, nr, nc = _bitpack_matrix(H)
    
    n_threads = nthreads()
    profiles = [fill(typemax(Int), nr + 1) for _ in 1:n_threads]
    thread_syndromes = [zeros(UInt64, max_wt + 1, num_words) for _ in 1:n_threads]

    @threads for c1 in 1:nc
        tid = threadid()
        syndromes = thread_syndromes[tid]
        profile = profiles[tid]
        
        syn_wt = 0
        @inbounds for w in 1:num_words
            val = H_packed[c1, w]
            syndromes[2, w] = val
            syn_wt += count_ones(val)
        end
        
        if syn_wt <= max_syndrome_wt
            @inbounds if 1 < profile[syn_wt + 1]
                profile[syn_wt + 1] = 1
            end
        end
        
        _confinement_dfs!(c1 + 1, 1, max_wt, nc, num_words, Δ, syndromes, H_packed, profile, max_syndrome_wt)
    end
    
    global_profile = fill(typemax(Int), nr + 1)
    for tid in 1:n_threads
        @inbounds for i in 1:nr+1
            if profiles[tid][i] < global_profile[i]
                global_profile[i] = profiles[tid][i]
            end
        end
    end
    
    conf_dict = Dict{Int, Int}()
    for (syn_wt_plus_1, min_e_wt) in enumerate(global_profile)
        if min_e_wt != typemax(Int) && syn_wt_plus_1 - 1 > 0
            conf_dict[syn_wt_plus_1 - 1] = min_e_wt
        end
    end
    
    return conf_dict
end

# -----------------------------------------------------------------------------
# 3. Verify Soundness
# -----------------------------------------------------------------------------
function _verify_QLTC_soundness_dfs!(start_col::Int, depth::Int, max_wt::Int, nc::Int, num_words::Int, 
                                Δ::Int, syndromes::Matrix{UInt64}, H_packed::Matrix{UInt64}, 
                                threshold::Float64, is_valid::Atomic{Int})
    if is_valid[] == 0 || depth == max_wt
        return
    end

    syn_wt = 0
    @inbounds for w in 1:num_words
        syn_wt += count_ones(syndromes[depth + 1, w])
    end

    best_future_ratio = Inf
    for k in (depth + 1):max_wt
        best_possible_syn_wt = max(0, syn_wt - (k - depth) * Δ)
        ratio_bound = best_possible_syn_wt / k
        if ratio_bound < best_future_ratio
            best_future_ratio = ratio_bound
        end
    end

    @inbounds if best_future_ratio >= threshold
        return
    end

    for c in start_col:nc
        next_syn_wt = 0
        @inbounds for w in 1:num_words
            val = syndromes[depth + 1, w] ⊻ H_packed[c, w]
            syndromes[depth + 2, w] = val
            next_syn_wt += count_ones(val)
        end

        if next_syn_wt > 0
            ratio = next_syn_wt / (depth + 1)
            @inbounds if ratio < threshold
                atomic_xchg!(is_valid, 0)
                return
            end
        end

        _verify_QLTC_soundness_dfs!(c + 1, depth + 1, max_wt, nc, num_words, Δ, syndromes, H_packed, threshold, is_valid)
        
        if is_valid[] == 0
            return
        end
    end
end

"""
    verify_QLTC_soundness(S::AbstractSubsystemCode, threshold::Float64, check_type::Symbol=:X; max_wt::Int=4)

Returns `true` if the code's soundness is strictly `>= threshold` for all errors up to `max_wt`.
Returns `false` and instantly aborts all threads the moment a violation is found.
"""
function verify_QLTC_soundness(S::AbstractSubsystemCode, threshold::Float64, check_type::Symbol=:X; max_wt::Int=4)
    check_type ∈ (:X, :Z) || throw(ArgumentError("check_type must be :X or :Z"))
    H_orig = check_type == :X ? X_stabilizers(S) : Z_stabilizers(S)
    
    col_wts = [count(!iszero, H_orig[:, c]) for c in 1:size(H_orig, 2)]
    p = sortperm(col_wts)
    H = H_orig[:, p]
    Δ = maximum(col_wts)
    
    H_packed, num_words, nr, nc = _bitpack_matrix(H)
    
    n_threads = nthreads()
    is_valid = Atomic{Int}(1)
    thread_syndromes = [zeros(UInt64, max_wt + 1, num_words) for _ in 1:n_threads]

    @threads for c1 in 1:nc
        if is_valid[] == 0
            continue
        end
        
        tid = threadid()
        syndromes = thread_syndromes[tid]
        
        syn_wt = 0
        @inbounds for w in 1:num_words
            val = H_packed[c1, w]
            syndromes[2, w] = val
            syn_wt += count_ones(val)
        end
        
        if syn_wt > 0
            ratio = syn_wt / 1.0
            if ratio < threshold
                atomic_xchg!(is_valid, 0)
                continue
            end
        end
        
        _verify_QLTC_soundness_dfs!(c1 + 1, 1, max_wt, nc, num_words, Δ, syndromes, H_packed, threshold, is_valid)
    end
    
    return is_valid[] == 1
end

# -----------------------------------------------------------------------------
# 4. Verify Confinement
# -----------------------------------------------------------------------------
function _verify_confinement_dfs!(start_col::Int, depth::Int, max_wt::Int, nc::Int, num_words::Int, 
                                  Δ::Int, syndromes::Matrix{UInt64}, H_packed::Matrix{UInt64}, 
                                  target_syndrome_wt::Int, is_confined::Atomic{Int})
    if is_confined[] == 0 || depth == max_wt
        return
    end

    syn_wt = 0
    @inbounds for w in 1:num_words
        syn_wt += count_ones(syndromes[depth + 1, w])
    end

    best_possible_future_syn_wt = max(0, syn_wt - (max_wt - depth) * Δ)
    if best_possible_future_syn_wt > target_syndrome_wt
        return
    end

    for c in start_col:nc
        next_syn_wt = 0
        @inbounds for w in 1:num_words
            val = syndromes[depth + 1, w] ⊻ H_packed[c, w]
            syndromes[depth + 2, w] = val
            next_syn_wt += count_ones(val)
        end

        if next_syn_wt <= target_syndrome_wt
            atomic_xchg!(is_confined, 0)
            return
        end

        _verify_confinement_dfs!(c + 1, depth + 1, max_wt, nc, num_words, Δ, syndromes, H_packed, target_syndrome_wt, is_confined)
        
        if is_confined[] == 0
            return
        end
    end
end

"""
    verify_confinement(S::AbstractSubsystemCode, target_syndrome_wt::Int, check_type::Symbol=:X; max_wt::Int=4)

Returns `true` if every error up to `max_wt` produces a syndrome weight strictly `> target_syndrome_wt`.
Returns `false` and instantly aborts if an error is found that fails to trigger the required syndrome weight.
"""
function verify_confinement(S::AbstractSubsystemCode, target_syndrome_wt::Int, check_type::Symbol=:X; max_wt::Int=4)
    check_type ∈ (:X, :Z) || throw(ArgumentError("check_type must be :X or :Z"))
    H_orig = check_type == :X ? X_stabilizers(S) : Z_stabilizers(S)
    
    col_wts = [count(!iszero, H_orig[:, c]) for c in 1:size(H_orig, 2)]
    p = sortperm(col_wts)
    H = H_orig[:, p]
    Δ = maximum(col_wts)
    
    H_packed, num_words, nr, nc = _bitpack_matrix(H)
    
    n_threads = nthreads()
    is_confined = Atomic{Int}(1)
    thread_syndromes = [zeros(UInt64, max_wt + 1, num_words) for _ in 1:n_threads]

    @threads for c1 in 1:nc
        if is_confined[] == 0
            continue
        end
        
        tid = threadid()
        syndromes = thread_syndromes[tid]
        
        syn_wt = 0
        @inbounds for w in 1:num_words
            val = H_packed[c1, w]
            syndromes[2, w] = val
            syn_wt += count_ones(val)
        end
        
        if syn_wt <= target_syndrome_wt
            atomic_xchg!(is_confined, 0)
            continue
        end
        
        _verify_confinement_dfs!(c1 + 1, 1, max_wt, nc, num_words, Δ, syndromes, H_packed, target_syndrome_wt, is_confined)
    end
    
    return is_confined[] == 1
end

"""
    normalized_laplacian_matrix(G::Grphs.SimpleGraph)

Return the normalized Laplacian matrix of the graph.
Essential for evaluating the expansion of irregular Tanner graphs.
"""
function normalized_laplacian_matrix(G::Grphs.SimpleGraph)
    A = Grphs.adjacency_matrix(G)
    d = Grphs.degree(G)
    
    # Safely handle isolated vertices to avoid division by zero
    inv_sqrt_d = [v > 0 ? 1.0 / sqrt(v) : 0.0 for v in d]
    D_inv_sqrt = spdiagm(0 => inv_sqrt_d)
    I_mat = spdiagm(0 => ones(Grphs.nv(G)))
    
    return I_mat - D_inv_sqrt * A * D_inv_sqrt
end

"""
    normalized_spectral_gap(G::Grphs.SimpleGraph)

Return the algebraic connectivity (second smallest eigenvalue) of the 
normalized Laplacian matrix. Bounded between 0 and 2.
"""
function normalized_spectral_gap(G::Grphs.SimpleGraph)
    L_norm = normalized_laplacian_matrix(G)
    evals = eigvals(Symmetric(Matrix(L_norm)))
    return length(evals) > 1 ? evals[2] : 0.0
end

"""
    fiedler_vector(S::AbstractSubsystemCode, check_type::Symbol=:both)

Returns the Fiedler vector (the eigenvector corresponding to the second smallest 
Laplacian eigenvalue). Useful for spectral graph partitioning and hardware layout.
Results are cached.
"""
function fiedler_vector(S::AbstractSubsystemCode, check_type::Symbol=:both)
    if CSSTrait(typeof(S)) == IsNotCSS()
        check_type == :both || @warn "check_type ignored for non-CSS codes. Using full graph."
        haskey(S.cache, :fiedler_vector) && return S.cache[:fiedler_vector]
        G = Tanner_graph(S)[1]
    else
        check_type ∈ (:X, :Z, :both) || throw(ArgumentError("check_type must be :X, :Z, or :both"))
        cache_key = Symbol("fiedler_vector_", check_type)
        haskey(S.cache, cache_key) && return S.cache[cache_key]
        
        G = check_type == :X ? Tanner_graph_X(S)[1] : (check_type == :Z ? Tanner_graph_Z(S)[1] : Tanner_graph(S)[1])
    end
    
    L = laplacian_matrix(G)
    evals, evecs = eigen(Symmetric(Matrix(L)))
    
    # Eigen returns sorted eigenvalues. Index 2 is the algebraic connectivity.
    vec = size(evecs, 2) > 1 ? evecs[:, 2] : zeros(Grphs.nv(G))
    
    if CSSTrait(typeof(S)) == IsNotCSS()
        S.cache[:fiedler_vector] = vec
    else
        S.cache[cache_key] = vec
    end
    return vec
end

"""
    is_topologically_connected(S::AbstractSubsystemCode, check_type::Symbol=:both)

Returns `true` if the underlying Tanner graph consists of a single connected component.
"""
function is_topologically_connected(S::AbstractSubsystemCode, check_type::Symbol=:both)
    if CSSTrait(typeof(S)) == IsNotCSS()
        G = Tanner_graph(S)[1]
    else
        check_type ∈ (:X, :Z, :both) || throw(ArgumentError("check_type must be :X, :Z, or :both"))
        G = check_type == :X ? Tanner_graph_X(S)[1] : (check_type == :Z ? Tanner_graph_Z(S)[1] : Tanner_graph(S)[1])
    end
    
    # We can use the native Graphs.jl function for maximum speed
    return Grphs.is_connected(G)
end

"""
    tanner_distance_bound(S::AbstractSubsystemCode, check_type::Symbol=:both)

Returns the spectral lower bound on the minimum distance of the code based on the 
second-largest adjacency eigenvalue of its Tanner graph.
Throws an error if the specified graph is not strictly regular.
"""
function tanner_distance_bound(S::AbstractSubsystemCode, check_type::Symbol=:both)
    if CSSTrait(typeof(S)) == IsNotCSS()
        G = Tanner_graph(S)[1]
        n_qubits = length(S)
    else
        check_type ∈ (:X, :Z, :both) || throw(ArgumentError("check_type must be :X, :Z, or :both"))
        G = check_type == :X ? Tanner_graph_X(S)[1] : (check_type == :Z ? Tanner_graph_Z(S)[1] : Tanner_graph(S)[1])
        n_qubits = length(S)
    end
    
    Grphs.is_regular(G) || throw(ArgumentError("The Tanner distance bound is only mathematically valid for regular graphs."))
    
    Δ = Grphs.degree(G)[1]
    
    A = Grphs.adjacency_matrix(G)
    evals = eigvals(Symmetric(Matrix(A)))
    
    # Get the second largest eigenvalue (λ)
    max_eval = maximum(abs.(evals))
    non_trivial_evals = filter(x -> abs(x) < max_eval - 1e-7, evals)
    λ = isempty(non_trivial_evals) ? 0.0 : maximum(abs.(non_trivial_evals))
    
    # Tanner bound equation
    bound = n_qubits * (Δ - λ) / (Δ^2 - 2λ + Δ)
    return max(1, floor(Int, bound))
end

"""
    estimated_edge_expansion(G::Grphs.SimpleGraph)

Estimates the edge expansion (isoperimetric number) of the graph using a Fiedler vector sweep.
Returns an upper bound on the true edge expansion.
"""
function estimated_edge_expansion(G::Grphs.SimpleGraph)
    L = laplacian_matrix(G)
    evals, evecs = eigen(Symmetric(Matrix(L)))
    v = size(evecs, 2) > 1 ? evecs[:, 2] : zeros(Grphs.nv(G))
    p = sortperm(v)
    
    min_exp = Inf
    n = Grphs.nv(G)
    
    S_set = falses(n)
    cut_size = 0
    
    for i in 1:(n ÷ 2)
        node = p[i]
        S_set[node] = true
        
        # Fast update of the cut size
        for neighbor in Grphs.neighbors(G, node)
            if S_set[neighbor]
                cut_size -= 1
            else
                cut_size += 1
            end
        end
        
        exp = cut_size / i
        if exp < min_exp
            min_exp = exp
        end
    end
    
    return min_exp == Inf ? 0.0 : min_exp
end

"""
    estimated_vertex_expansion(G::Grphs.SimpleGraph)

Estimates the vertex expansion of the graph using a Fiedler vector sweep.
Returns an upper bound on the true vertex expansion.
"""
function estimated_vertex_expansion(G::Grphs.SimpleGraph)
    L = laplacian_matrix(G)
    evals, evecs = eigen(Symmetric(Matrix(L)))
    v = size(evecs, 2) > 1 ? evecs[:, 2] : zeros(Grphs.nv(G))
    p = sortperm(v)
    
    min_exp = Inf
    n = Grphs.nv(G)
    
    S_set = falses(n)
    in_boundary = falses(n)
    boundary_size = 0
    
    for i in 1:(n ÷ 2)
        node = p[i]
        S_set[node] = true
        if in_boundary[node]
            boundary_size -= 1 # Node moved from boundary to interior
        end
        
        for neighbor in Grphs.neighbors(G, node)
            if !S_set[neighbor] && !in_boundary[neighbor]
                in_boundary[neighbor] = true
                boundary_size += 1
            end
        end
        
        exp = boundary_size / i
        if exp < min_exp
            min_exp = exp
        end
    end
    
    return min_exp == Inf ? 0.0 : min_exp
end

"""
    edge_expansion_bounds(G::Grphs.SimpleGraph)

Returns a tuple `(lower_bound, upper_bound)` for the edge expansion of the graph, 
computed strictly using the Cheeger inequalities on the Laplacian eigenvalues.
"""
function edge_expansion_bounds(G::Grphs.SimpleGraph)
    L = laplacian_matrix(G)
    evals = eigvals(Symmetric(Matrix(L)))
    
    # Algebraic connectivity (second smallest Laplacian eigenvalue)
    μ_2 = length(evals) > 1 ? evals[2] : 0.0
    
    # Max degree
    Δ = maximum(Grphs.degree(G))
    
    lower_bound = μ_2 / 2.0
    upper_bound = sqrt(2 * Δ * μ_2)
    
    return (lower_bound, upper_bound)
end

"""
    estimated_bipartite_vertex_expansion(S::AbstractSubsystemCode, check_type::Symbol=:both; max_subset_fraction::Float64=0.5)

Estimates the bipartite vertex expansion (qubits to checks) of the specified sub-graph 
(`:X`, `:Z`, or `:both` for non-CSS) of the code `S`.
`max_subset_fraction` controls the maximum size of the qubit subset `S` as a fraction of total qubits.
Results are cached dynamically based on the fraction provided.
"""
function estimated_bipartite_vertex_expansion(S::AbstractSubsystemCode, check_type::Symbol=:both; max_subset_fraction::Float64=0.5)
    if CSSTrait(typeof(S)) == IsNotCSS()
        check_type == :both || @warn "check_type ignored for non-CSS codes. Using full graph."
        
        # Dynamic cache key combining function name and fraction
        cache_key = Symbol("estimated_bipartite_vertex_expansion_", max_subset_fraction)
        haskey(S.cache, cache_key) && return S.cache[cache_key]
        
        G_tuple = Tanner_graph(S)
    else
        check_type ∈ (:X, :Z, :both) || throw(ArgumentError("check_type must be :X, :Z, or :both"))
        
        cache_key = Symbol("estimated_bipartite_vertex_expansion_", check_type, "_", max_subset_fraction)
        haskey(S.cache, cache_key) && return S.cache[cache_key]
        
        G_tuple = check_type == :X ? Tanner_graph_X(S) : (check_type == :Z ? Tanner_graph_Z(S) : Tanner_graph(S))
    end
    
    G = G_tuple[1]
    n_qubits = length(G_tuple[2]) 
    max_subset_size = floor(Int, n_qubits * max_subset_fraction)
    
    val = estimated_bipartite_vertex_expansion(G, n_qubits, max_subset_size=max_subset_size)
    S.cache[cache_key] = val
    return val
end

for func in (:estimated_edge_expansion, :estimated_vertex_expansion, :edge_expansion_bounds)
    @eval begin
        """
            $($func)(S::AbstractSubsystemCode, check_type::Symbol=:both)
        
        Computes the $($func) of the specified sub-graph (`:X`, `:Z`, or `:both` for non-CSS) of the code `S`.
        Results are cached.
        """
        function $func(S::AbstractSubsystemCode, check_type::Symbol=:both)
            if CSSTrait(typeof(S)) == IsNotCSS()
                check_type == :both || @warn "check_type ignored for non-CSS codes. Using full graph."
                
                cache_key = Symbol($func)
                haskey(S.cache, cache_key) && return S.cache[cache_key]
                
                G = Tanner_graph(S)[1]
            else
                check_type ∈ (:X, :Z, :both) || throw(ArgumentError("check_type must be :X, :Z, or :both"))
                
                cache_key = Symbol($func, "_", check_type)
                haskey(S.cache, cache_key) && return S.cache[cache_key]
                
                G = check_type == :X ? Tanner_graph_X(S)[1] : (check_type == :Z ? Tanner_graph_Z(S)[1] : Tanner_graph(S)[1])
            end
            
            val = $func(G) # Dispatch to the base Graph function
            S.cache[cache_key] = val
            return val
        end
    end
end

using Base.Threads: @threads, nthreads, threadid, Atomic, atomic_xchg!

function _verify_expander_dfs!(start_col::Int, depth::Int, max_wt::Int, nc::Int, num_words::Int, 
                               neighborhoods::Matrix{UInt64}, H_packed::Matrix{UInt64}, 
                               A::Float64, is_valid::Atomic{Int})
    if is_valid[] == 0 || depth == max_wt
        return
    end

    for c in start_col:nc
        next_syn_wt = 0
        @inbounds for w in 1:num_words
            # Boolean OR (|) instead of XOR (⊻) to get true topological neighborhood
            val = neighborhoods[depth + 1, w] | H_packed[c, w]
            neighborhoods[depth + 2, w] = val
            next_syn_wt += count_ones(val)
        end

        # If it fails the expansion threshold, abort all threads
        if next_syn_wt < A * (depth + 1)
            atomic_xchg!(is_valid, 0)
            return
        end

        # TOPOLOGICAL PRUNING:
        # Since weight monotonically increases, if we have already amassed enough 
        # neighbors to satisfy the maximum possible future threshold (A * max_wt), 
        # no superset down this branch can possibly fail. Prune it!
        if next_syn_wt >= A * max_wt
            continue
        end

        _verify_expander_dfs!(c + 1, depth + 1, max_wt, nc, num_words, neighborhoods, H_packed, A, is_valid)
        
        if is_valid[] == 0
            return
        end
    end
end

"""
    is_expander(H::CTMatrixTypes, γ::Float64, A::Float64)

Returns `true` if the bipartite graph represented by parity-check matrix `H` 
is a `(γ, A)`-expander. This strictly verifies that every subset of columns `S` 
with size up to `γ * nc` has a topological neighborhood of at least `A * |S|`.
Uses a highly optimized, multithreaded DFS with Bitwise-OR logic.
"""
function is_expander(H::CTMatrixTypes, γ::Float64, A::Float64)
    nr, nc = size(H)
    max_wt = floor(Int, γ * nc)
    if max_wt < 1
        return true
    end
    
    # SORTING HEURISTIC: Sort columns descending by weight so the topological 
    # neighborhood grows as fast as possible, triggering the B&B cutoff early.
    col_wts = [count(!iszero, H[:, c]) for c in 1:nc]
    p = sortperm(col_wts, rev=true)
    H_sorted = H[:, p]
    
    H_packed, num_words, _, _ = _bitpack_matrix(H_sorted)
    
    n_threads = nthreads()
    is_valid = Atomic{Int}(1)
    thread_neighborhoods = [zeros(UInt64, max_wt + 1, num_words) for _ in 1:n_threads]

    @threads for c1 in 1:nc
        if is_valid[] == 0
            continue
        end
        
        tid = threadid()
        neighborhoods = thread_neighborhoods[tid]
        
        syn_wt = 0
        @inbounds for w in 1:num_words
            val = H_packed[c1, w]
            neighborhoods[2, w] = val
            syn_wt += count_ones(val)
        end
        
        if syn_wt < A * 1.0
            atomic_xchg!(is_valid, 0)
            continue
        end
        
        if syn_wt >= A * max_wt
            continue # Prune instantly!
        end
        
        _verify_expander_dfs!(c1 + 1, 1, max_wt, nc, num_words, neighborhoods, H_packed, A, is_valid)
    end
    
    return is_valid[] == 1
end

"""
    is_expander(S::AbstractSubsystemCode, γ::Float64, A::Float64, check_type::Symbol=:both)

Returns `true` if the specified sub-graph (`:X`, `:Z`, or `:both` for non-CSS) 
of the code `S` acts as a `(γ, A)`-expander from qubits to checks.
"""
function is_expander(S::AbstractSubsystemCode, γ::Float64, A::Float64, check_type::Symbol=:both)
    if CSSTrait(typeof(S)) == IsNotCSS()
        check_type == :both || @warn "check_type ignored for non-CSS codes. Using full matrix."
        H = stabilizers(S)
    else
        check_type ∈ (:X, :Z, :both) || throw(ArgumentError("check_type must be :X, :Z, or :both"))
        if check_type == :X
            H = X_stabilizers(S)
        elseif check_type == :Z
            H = Z_stabilizers(S)
        else
            # For the full tripartite graph, the neighborhood is all X and Z checks combined
            H = vcat(X_stabilizers(S), Z_stabilizers(S))
        end
    end
    
    return is_expander(H, γ, A)
end

"""
    is_left_right_expander(H::CTMatrixTypes, γ_L::Float64, A_L::Float64, γ_R::Float64, A_R::Float64)

Verifies if the bipartite graph represented by biadjacency matrix `H` is a two-sided left-right expander.
Checks if columns (Left) expand into rows (Right) with parameters `(γ_L, A_L)`, AND 
if rows (Right) expand into columns (Left) with parameters `(γ_R, A_R)`.
"""
function is_left_right_expander(H::CTMatrixTypes, γ_L::Float64, A_L::Float64, γ_R::Float64, A_R::Float64)
    # 1. Check Left-to-Right expansion
    is_expander(H, γ_L, A_L) || return false
    
    # 2. Extract transposed matrix for Right-to-Left check.
    # We materialize a dense Int matrix to guarantee our bit-packer 
    # runs at maximum speed regardless of the underlying CTMatrixTypes format.
    nr, nc = size(H)
    H_T = zeros(Int, nc, nr)
    for r in 1:nr
        for c in 1:nc
            if !iszero(H[r, c])
                H_T[c, r] = 1
            end
        end
    end
    
    # Check Right-to-Left expansion
    return is_expander(H_T, γ_R, A_R)
end

"""
    is_bipartite_expander(G::Grphs.SimpleGraph, left::Vector{Int}, right::Vector{Int}, γ::Float64, A::Float64)

Returns `true` if the bipartition defined by `left` and `right` nodes in `G` forms a one-sided 
`(γ, A)`-expander strictly from the `left` set to the `right` set.
"""
function is_bipartite_expander(G::Grphs.SimpleGraph, left::Vector{Int}, right::Vector{Int}, γ::Float64, A::Float64)
    nr = length(right)
    nc = length(left)
    H = zeros(Int, nr, nc)
    
    # Map right nodes to row indices for O(1) lookup
    right_map = Dict(v => i for (i, v) in enumerate(right))
    
    for (c, u) in enumerate(left)
        for v in Grphs.neighbors(G, u)
            if haskey(right_map, v)
                H[right_map[v], c] = 1
            end
        end
    end
    
    return is_expander(H, γ, A)
end

"""
    is_left_right_expander(G::Grphs.SimpleGraph, left::Vector{Int}, right::Vector{Int}, 
                           γ_L::Float64, A_L::Float64, γ_R::Float64, A_R::Float64)

Verifies if the bipartition in `G` forms a two-sided left-right expander.
"""
function is_left_right_expander(G::Grphs.SimpleGraph, left::Vector{Int}, right::Vector{Int}, 
                                γ_L::Float64, A_L::Float64, γ_R::Float64, A_R::Float64)
    # Check Left -> Right
    is_bipartite_expander(G, left, right, γ_L, A_L) || return false
    
    # Check Right -> Left
    return is_bipartite_expander(G, right, left, γ_R, A_R)
end

"""
    is_left_right_expander(C::AbstractLinearCode, γ_L::Float64, A_L::Float64, γ_R::Float64, A_R::Float64)

Verifies if the classical Tanner graph of code `C` is a two-sided left-right expander, 
where bits (Left) expand to checks (Right), and checks (Right) expand to bits (Left).
"""
function is_left_right_expander(C::AbstractLinearCode, γ_L::Float64, A_L::Float64, γ_R::Float64, A_R::Float64)
    H = parity_check_matrix(C)
    return is_left_right_expander(H, γ_L, A_L, γ_R, A_R)
end

"""
    verify_quantum_tanner_structure(S::AbstractSubsystemCode, 
                                    G_A::Grphs.SimpleGraph, G_B::Grphs.SimpleGraph, 
                                    C_A::AbstractLinearCode, C_B::AbstractLinearCode)

Verifies the mathematical requirements for a valid Quantum Tanner Code construction.
Returns `true` if all structural, topological, and spectral conditions are met.

# Checks:
1. The global code is a valid CSS code.
2. The base graphs `G_A` and `G_B` are regular.
3. The degrees of the base graphs exactly match the block lengths of the local classical codes.
4. The Sipser-Spielman spectral condition: The classical minimum distances `d_A` and `d_B` 
   strictly exceed the second-largest adjacency eigenvalues `λ_A` and `λ_B` of their respective graphs.
"""
function verify_quantum_tanner_structure(S::AbstractSubsystemCode, 
                                         G_A::Grphs.SimpleGraph, G_B::Grphs.SimpleGraph, 
                                         C_A::AbstractLinearCode, C_B::AbstractLinearCode)
    
    is_valid = true

    # 1. CSS Verification
    if CSSTrait(typeof(S)) == IsNotCSS()
        @warn "Quantum Tanner codes must be CSS codes."
        is_valid = false
    end

    # 2. Graph Regularity Check
    if !Grphs.is_regular(G_A) || !Grphs.is_regular(G_B)
        @warn "Base graphs G_A and G_B must be regular."
        is_valid = false
    end

    # Fast abort if non-regular, because degree and eigenvalue metrics will fail or be meaningless
    if !is_valid
        return false
    end

    Δ_A = Grphs.degree(G_A)[1]
    Δ_B = Grphs.degree(G_B)[1]

    # 3. Block Length Matching
    if Δ_A != length(C_A)
        @warn "Degree of G_A ($Δ_A) does not match the block length of local code C_A ($(length(C_A)))."
        is_valid = false
    end
    if Δ_B != length(C_B)
        @warn "Degree of G_B ($Δ_B) does not match the block length of local code C_B ($(length(C_B)))."
        is_valid = false
    end

    # 4. The Spectral/Distance Bound (d_0 > λ)
    # Note: Our spectral_gap function returns the second-largest adjacency eigenvalue (λ).
    λ_A = spectral_gap(G_A)
    λ_B = spectral_gap(G_B)
    
    d_A = minimum_distance(C_A)
    d_B = minimum_distance(C_B)

    if d_A <= λ_A
        @warn "Spectral bound failed for A-complex: Local distance d_A ($d_A) must be strictly greater than λ_A ($λ_A)."
        is_valid = false
    else
        @info "A-complex spectral bound passed: d_A ($d_A) > λ_A ($λ_A)"
    end

    if d_B <= λ_B
        @warn "Spectral bound failed for B-complex: Local distance d_B ($d_B) must be strictly greater than λ_B ($λ_B)."
        is_valid = false
    else
        @info "B-complex spectral bound passed: d_B ($d_B) > λ_B ($λ_B)"
    end

    if is_valid
        @info "Quantum Tanner structure verified successfully."
    end

    return is_valid
end

"""
    verify_quantum_tanner_structure(S::AbstractSubsystemCode, G_base::Grphs.SimpleGraph, C_local::AbstractLinearCode)

Symmetric convenience wrapper for Quantum Tanner Codes constructed from a single base graph and a single local classical code.
"""
function verify_quantum_tanner_structure(S::AbstractSubsystemCode, G_base::Grphs.SimpleGraph, C_local::AbstractLinearCode)
    return verify_quantum_tanner_structure(S, G_base, G_base, C_local, C_local)
end

"""
    sipser_spielman_guarantees(H::CTMatrixTypes, γ::Float64, A::Float64)

Evaluates the Sipser-Spielman guarantees for a given parity-check matrix `H` 
that is known to be a `(γ, A)`-expander. 

Returns a NamedTuple containing the expansion factor `ϵ`, a boolean indicating 
if linear distance is guaranteed, the guaranteed minimum distance bound, and the 
guaranteed parallel bit-flip decoding radius.
"""
function sipser_spielman_guarantees(H::CTMatrixTypes, γ::Float64, A::Float64)
    nr, nc = size(H)
    
    # Check if the matrix is left-regular (all columns have the same weight)
    col_wts = [count(!iszero, H[:, c]) for c in 1:nc]
    Δ = col_wts[1]
    
    if any(w != Δ for w in col_wts)
        @warn "Sipser-Spielman bounds technically assume a strictly left-regular bipartite graph. Using average degree."
        Δ = sum(col_wts) / nc
    end
    
    # ϵ is the fraction of the maximum possible neighborhood
    ϵ = A / Δ
    
    # Maximum size of sets S that achieve this expansion
    max_S = floor(Int, γ * nc)
    
    # Guarantee 1: Unique Neighbors & Distance (ϵ > 1/2)
    has_linear_distance = ϵ > 0.5
    guaranteed_d = has_linear_distance ? max_S + 1 : 0
    
    # Guarantee 2: Parallel Bit-Flip Decoding (ϵ > 3/4)
    has_guaranteed_decoding = ϵ > 0.75
    bit_flip_radius = has_guaranteed_decoding ? floor(Int, max_S / 2) : 0
    
    return (
        epsilon = ϵ,
        delta = Δ,
        has_linear_distance = has_linear_distance,
        guaranteed_distance = guaranteed_d,
        has_guaranteed_decoding = has_guaranteed_decoding,
        bit_flip_radius = bit_flip_radius
    )
end

"""
    sipser_spielman_guarantees(S::AbstractSubsystemCode, γ::Float64, A::Float64, check_type::Symbol=:both)

Evaluates the Sipser-Spielman guarantees for the specified sub-graph of a quantum code.
"""
function sipser_spielman_guarantees(S::AbstractSubsystemCode, γ::Float64, A::Float64, check_type::Symbol=:both)
    if CSSTrait(typeof(S)) == IsNotCSS()
        H = stabilizers(S)
    else
        check_type ∈ (:X, :Z, :both) || throw(ArgumentError("check_type must be :X, :Z, or :both"))
        if check_type == :X
            H = X_stabilizers(S)
        elseif check_type == :Z
            H = Z_stabilizers(S)
        else
            H = vcat(X_stabilizers(S), Z_stabilizers(S))
        end
    end
    
    return sipser_spielman_guarantees(H, γ, A)
end

# # Example workflow the user would run:
# # 1. Verify the code expands at 0.05 (5% of qubits) with an expansion multiplier of 3.8
# if is_expander(H, 0.05, 3.8)
    
#     # 2. Extract the physical bounds
#     bounds = sipser_spielman_guarantees(H, 0.05, 3.8)
#     println("Expansion Factor ϵ: ", bounds.epsilon)
#     println("Guaranteed Distance: ", bounds.guaranteed_distance)
# end

"""
    evaluate_single_shot_soundness(H::CTMatrixTypes, t::Int; max_error_wt::Int=5)

Evaluates the Campbell (QEC) single-shot soundness of the parity-check matrix `H`.
Returns a dictionary mapping syndrome weight `w` (for `w < t`) to the maximum minimum-weight 
error `E^{red}` required to produce it. This directly defines the bounding function `f(w)`.

Uses a highly optimized Breadth-First Search (BFS) over bit-packed matrices.
"""
function evaluate_single_shot_soundness(H::CTMatrixTypes, t::Int; max_error_wt::Int=5)
    nc = size(H, 2)
    H_packed, num_words, _, _ = _bitpack_matrix(H)
    
    # Track visited syndromes. In BFS, first visit == absolute minimum error weight.
    visited = Set{Vector{UInt64}}()
    
    # Initialize the BFS queue. We use double-buffering (current and next level arrays) 
    # which is significantly faster in Julia than pushing/popping from a standard Queue.
    # State tuple: (syndrome_array, last_column_flipped)
    current_level = [(zeros(UInt64, num_words), 0)]
    push!(visited, current_level[1][1])
    
    # f_vals maps syndrome weight -> max E^{red} weight
    f_vals = Dict{Int, Int}(0 => 0) # Strictly enforced: f(0) = 0
    
    for err_wt in 1:max_error_wt
        next_level = Tuple{Vector{UInt64}, Int}[]
        
        for (syn, last_col) in current_level
            # Only flip columns > last_col to perfectly avoid permutation redundancies 
            # (e.g., flipping col 1 then 2 is identical to 2 then 1)
            for c in (last_col + 1):nc
                new_syn = copy(syn)
                syn_wt = 0
                
                # Bitwise XOR to apply the error
                @inbounds for w in 1:num_words
                    new_syn[w] ⊻= H_packed[c, w]
                    syn_wt += count_ones(new_syn[w])
                end
                
                # If we haven't seen this syndrome yet, this err_wt is its absolute E^{red}
                if !(new_syn in visited)
                    push!(visited, new_syn)
                    push!(next_level, (new_syn, c))
                    
                    # If it's within our single-shot evaluation boundary, update the bounding function
                    if syn_wt < t
                        f_vals[syn_wt] = max(get(f_vals, syn_wt, 0), err_wt)
                    end
                end
            end
        end
        
        current_level = next_level
        # Early exit if we've exhausted the reachable syndrome space
        if isempty(current_level)
            break
        end
    end
    
    return f_vals
end

"""
    evaluate_single_shot_soundness(S::AbstractSubsystemCode, t::Int, check_type::Symbol=:both; max_error_wt::Int=5)

Evaluates the Campbell (QEC) single-shot soundness for the specified stabilizers of a quantum code.
Returns the bounding function `f(w)` as a dictionary mapping syndrome weights (up to `t-1`) 
to their worst-case minimum error weight.
"""
function evaluate_single_shot_soundness(S::AbstractSubsystemCode, t::Int, check_type::Symbol=:both; max_error_wt::Int=5)
    if CSSTrait(typeof(S)) == IsNotCSS()
        H = stabilizers(S)
    else
        check_type ∈ (:X, :Z, :both) || throw(ArgumentError("check_type must be :X, :Z, or :both"))
        if check_type == :X
            H = X_stabilizers(S)
        elseif check_type == :Z
            H = Z_stabilizers(S)
        else
            H = vcat(X_stabilizers(S), Z_stabilizers(S))
        end
    end
    
    return evaluate_single_shot_soundness(H, t, max_error_wt=max_error_wt)
end

"""
    evaluate_confinement(H::CTMatrixTypes, t::Int)

Evaluates the (t, f)-confinement of the parity-check matrix `H`.
Returns a dictionary representing the function `f(s)`, which maps a syndrome weight `s` 
to the maximum error weight `w` (where `w <= t`) that produced it.

A code with good confinement partitions low-energy states into well-separated clusters, 
keeping residual errors from minimum-weight decoders bounded.
"""
function evaluate_confinement(H::CTMatrixTypes, t::Int)
    nc = size(H, 2)
    H_packed, num_words, _, _ = _bitpack_matrix(H)
    
    # f_vals maps syndrome weight -> maximum error weight that caused it
    f_vals = Dict{Int, Int}(0 => 0) # f(0) = 0 is strictly required
    
    # BFS Queue state: (syndrome_array, last_col_flipped)
    current_level = [(zeros(UInt64, num_words), 0)]
    
    for err_wt in 1:t
        next_level = Tuple{Vector{UInt64}, Int}[]
        
        for (syn, last_col) in current_level
            # Only flip columns strictly greater than last_col to avoid permutations
            for c in (last_col + 1):nc
                new_syn = copy(syn)
                syn_wt = 0
                
                @inbounds for w in 1:num_words
                    new_syn[w] ⊻= H_packed[c, w]
                    syn_wt += count_ones(new_syn[w])
                end
                
                # Update the bounding function: f(syn_wt) >= err_wt
                if !haskey(f_vals, syn_wt)
                    f_vals[syn_wt] = err_wt
                else
                    f_vals[syn_wt] = max(f_vals[syn_wt], err_wt)
                end
                
                push!(next_level, (new_syn, c))
            end
        end
        
        current_level = next_level
        if isempty(current_level)
            break
        end
    end
    
    return f_vals
end

"""
    evaluate_confinement(S::AbstractSubsystemCode, t::Int, check_type::Symbol=:both)

Evaluates the (t, f)-confinement bounds for the specified stabilizers of a quantum code `S`.
"""
function evaluate_confinement(S::AbstractSubsystemCode, t::Int, check_type::Symbol=:both)
    if CSSTrait(typeof(S)) == IsNotCSS()
        H = stabilizers(S)
    else
        check_type ∈ (:X, :Z, :both) || throw(ArgumentError("check_type must be :X, :Z, or :both"))
        if check_type == :X
            H = X_stabilizers(S)
        elseif check_type == :Z
            H = Z_stabilizers(S)
        else
            H = vcat(X_stabilizers(S), Z_stabilizers(S))
        end
    end
    
    return evaluate_confinement(H, t)
end

"""
    cosystolic_expansion(C::ChainComplex, k::Int; max_wt::Int=4, upper_bound::Float64=Inf)

Computes the exact cosystolic expansion ratio h_k for chains in degree `k` up to weight `max_wt`.
This relies on extracting the Oscar boundary maps and routing them to the optimized 
multithreaded Branch and Bound soundness engine.

Note: Currently strictly evaluates complexes over F_2.
"""
function cosystolic_expansion(C::ChainComplex, k::Int; max_wt::Int=4, upper_bound::Float64=Inf)
    # Check if the degree exists in our complex
    idx = findfirst(==(k), C.degrees)
    if idx === nothing
        throw(KeyError("Degree $k not found in the chain complex grading."))
    end
    
    ∂_k = C[k]
    
    # We also need ∂_{k+1} to establish the stabilizer/boundary equivalence class.
    # If k is the highest degree, the incoming map is 0, so the distance is just the Hamming weight.
    if idx == 1
        ∂_k_plus_1 = zero_matrix(base_ring(∂_k), ncols(∂_k), ncols(∂_k)) # Dummy empty matrix
    else
        deg_in = C.degrees[idx - 1]
        ∂_k_plus_1 = C[deg_in]
    end
    
    # 1. Convert Oscar matrices to standard dense Julia Int matrices for the bit-packer
    # (Assuming the base ring is F_2; if not, cosystolic expansion uses generalized Hamming weights)
    H_k = _oscar_to_dense_int(∂_k)
    
    # In our QEC soundness engine, the parity-check matrix (H) is the boundary map.
    # The solver intrinsically handles the topological expansion.
    # We bypass the `AbstractSubsystemCode` wrapper and call our DFS engine directly.
    
    nc = size(H_k, 2)
    
    # Sort columns descending by weight to trigger early Branch and Bound pruning
    col_wts = [count(!iszero, H_k[:, c]) for c in 1:nc]
    p = sortperm(col_wts, rev=true)
    H_sorted = H_k[:, p]
    Δ = maximum(col_wts)
    
    H_packed, num_words, nr, _ = _bitpack_matrix(H_sorted)
    
    # Launch the DFS soundness solver we wrote earlier
    n_threads = Base.Threads.nthreads()
    min_ratios = fill(upper_bound, n_threads)
    thread_syndromes = [zeros(UInt64, max_wt + 1, num_words) for _ in 1:n_threads]

    Base.Threads.@threads for c1 in 1:nc
        tid = Base.Threads.threadid()
        syndromes = thread_syndromes[tid]
        
        syn_wt = 0
        @inbounds for w in 1:num_words
            val = H_packed[c1, w]
            syndromes[2, w] = val
            syn_wt += count_ones(val)
        end
        
        if syn_wt > 0
            ratio = syn_wt / 1.0
            @inbounds if ratio < min_ratios[tid]
                min_ratios[tid] = ratio
            end
        end
        
        _soundness_dfs!(c1 + 1, 1, max_wt, nc, num_words, Δ, syndromes, H_packed, min_ratios, tid)
    end
    
    return minimum(min_ratios)
end

# Helper to bridge Oscar and our native bit-packer
function _oscar_to_dense_int(M)
    nr, nc = nrows(M), ncols(M)
    dense = zeros(Int, nr, nc)
    for r in 1:nr
        for c in 1:nc
            if !iszero(M[r, c])
                dense[r, c] = 1
            end
        end
    end
    return dense
end

"""
    verify_modular_expansion(H::Generic.MatSpaceElem, threshold::Float64; max_module_wt::Int=2)

Verifies the modular expansion of a parity-check matrix `H` defined over a group ring F_2[G].
Iterates directly through the module elements (vectors of group algebra elements) up to 
a block weight of `max_module_wt`.

Returns `true` if the ratio (syndrome block weight / error block weight) is >= `threshold`.
"""
function verify_modular_expansion(H::Generic.MatSpaceElem, threshold::Float64; max_module_wt::Int=2)
    R = base_ring(H) # This is the GroupAlgebra object in Oscar
    nr, nc = nrows(H), ncols(H)
    
    # We want to explore all non-zero module vectors of weight <= max_module_wt.
    # To do this natively, we need the underlying group elements to construct 
    # non-zero ring elements.
    G = group(R)
    group_elements = collect(G)
    
    # A module vector is an array of ring elements of length `nc`
    current_error = [zero(R) for _ in 1:nc]
    
    # Helper recursive worker for the algebraic search
    return _modular_expansion_dfs!(1, 0, max_module_wt, current_error, H, R, group_elements, nc, threshold)
end

function _modular_expansion_dfs!(start_idx::Int, current_wt::Int, max_wt::Int, 
                                 current_error::Vector{T}, H, R, group_elements, 
                                 nc::Int, threshold::Float64) where T
    
    # Evaluate the expansion ratio of the current non-zero module state
    if current_wt > 0
        # Compute the syndrome natively using Oscar's matrix-vector multiplication
        # We temporarily promote our vector to an Oscar column matrix
        err_matrix = matrix(R, nc, 1, current_error)
        syndrome = H * err_matrix
        
        # Calculate block weights (number of non-zero coordinates in the module vector)
        err_block_wt = current_wt
        syn_block_wt = 0
        for r in 1:nrows(syndrome)
            if !iszero(syndrome[r, 1])
                syn_block_wt += 1
            end
        end
        
        ratio = syn_block_wt / err_block_wt
        if ratio < threshold
            return false # Short-circuit immediately on failure
        end
    end
    
    if current_wt == max_wt
        return true
    end
    
    # Loop over module coordinates to inject errors
    for i in start_idx:nc
        # For a given coordinate, loop over non-zero group ring elements.
        # To keep it bounded, we look at single group element injections (monomials).
        for g in group_elements
            # Create the monomial element g in the group ring R
            g_elem = R(g)
            
            # Inject error into coordinate i
            current_error[i] = g_elem
            
            # Recurse to next weight level
            if !_modular_expansion_dfs!(i + 1, current_wt + 1, max_wt, current_error, H, R, group_elements, nc, threshold)
                return false
            end
            
            # Backtrack
            current_error[i] = zero(R)
        end
    end
    
    return true
end

import Oscar: base_ring, matrix, nrows, ncols, zero

"""
    verify_relative_modular_expansion(C::ChainComplex, threshold::Float64; max_module_wt::Int=2)

Evaluates the relative modular expansion of a quantum ChainComplex over F_2[G].
D1 (C.d[1]) maps errors to syndromes (H_X).
D2 (C.d[2]) maps stabilizer generators to errors (H_Z^T).

Returns `true` if min( syn_block_wt / dist_block(e, im(D2)) ) >= threshold.
"""
function verify_relative_modular_expansion(C::ChainComplex, threshold::Float64; max_module_wt::Int=2)
    D1 = C.d[1] # H_X: C_1 -> C_0
    D2 = C.d[2] # H_Z^T: C_2 -> C_1
    
    R = base_ring(D1)
    G = group(R)
    group_elements = collect(G)
    
    n_errors = ncols(D1) # Dimension of C_1
    
    # Pre-lift D2 to a flat binary matrix ONCE for fast distance calculations
    # (Assuming you have your `lift_matrix` utility available)
    D2_flat = lift_matrix(D2)
    
    current_error = [zero(R) for _ in 1:n_errors]
    
    return _relative_expansion_dfs!(
        1, 0, max_module_wt, current_error, 
        D1, D2_flat, R, group_elements, n_errors, threshold
    )
end

function _relative_expansion_dfs!(start_idx::Int, current_wt::Int, max_wt::Int, 
                                  current_error::Vector{T}, D1, D2_flat, R, group_elements, 
                                  nc::Int, threshold::Float64) where T
    
    if current_wt > 0
        err_matrix = matrix(R, nc, 1, current_error)
        syndrome = D1 * err_matrix
        
        # 1. Calculate Syndrome Block Weight
        syn_block_wt = 0
        for r in 1:nrows(syndrome)
            if !iszero(syndrome[r, 1])
                syn_block_wt += 1
            end
        end
        
        # 2. Calculate Relative Distance to Stabilizer Image
        # We only compute distance if the syndrome weight threatens our threshold.
        # If syn_block_wt is already large enough relative to absolute weight, 
        # it will definitely satisfy the threshold relative to the (smaller) distance.
        if (syn_block_wt / current_wt) < threshold
            
            # Map the group ring vector to a flat F_2 vector
            flat_error = _lift_vector(current_error, R)
            
            # Use classical min-dist solver on the coset: min wt(e + s) for s in im(D2_flat)
            # We then map the returned flat F_2 vector BACK to block weight.
            min_flat_err = find_minimum_coset_representative(flat_error, D2_flat)
            rel_dist = _calculate_block_weight(min_flat_err, length(group(R)))
            
            # Logical errors (rel_dist > 0, syn_wt == 0) will safely fail here
            if rel_dist > 0 && (syn_block_wt / rel_dist) < threshold
                return false
            end
        end
    end
    
    if current_wt == max_wt
        return true
    end
    
    for i in start_idx:nc
        for g in group_elements
            current_error[i] = R(g)
            
            if !_relative_expansion_dfs!(i + 1, current_wt + 1, max_wt, current_error, 
                                         D1, D2_flat, R, group_elements, nc, threshold)
                return false
            end
            
            current_error[i] = zero(R)
        end
    end
    
    return true
end

# --- Helpers to bridge Algebra and Geometry ---

"""
Converts a module vector over F_2[G] to a flat binary vector.
"""
function _lift_vector(e_vec::Vector{T}, R) where T
    G = group(R)
    N = length(G)
    flat_len = length(e_vec) * N
    flat = zeros(Int, flat_len)
    
    # Assuming group elements have a deterministic indexed order 1:N
    group_idx = Dict(g => i for (i, g) in enumerate(collect(G)))
    
    for (block_idx, poly) in enumerate(e_vec)
        if !iszero(poly)
            # Iterate through terms in the group ring polynomial
            for (coeff, g) in terms(poly) 
                if coeff == 1
                    flat[(block_idx - 1) * N + group_idx[g]] = 1
                end
            end
        end
    end
    return flat
end

"""
Calculates block weight directly from a flat binary vector.
"""
function _calculate_block_weight(flat_vec::Vector{Int}, lift_factor::Int)
    block_wt = 0
    num_blocks = length(flat_vec) ÷ lift_factor
    for b in 1:num_blocks
        start_idx = (b - 1) * lift_factor + 1
        end_idx = b * lift_factor
        if any(!iszero, @view flat_vec[start_idx:end_idx])
            block_wt += 1
        end
    end
    return block_wt
end

# import Oscar: base_ring, matrix, nrows, ncols, zero, group, terms

"""
    verify_relative_modular_expansion(C::ChainComplex, threshold::Float64; 
                                      max_module_wt::Int=2, exact::Bool=true, confidence::Float64=0.99)

Evaluates the relative modular expansion of a quantum ChainComplex over F_2[G].
D1 (C.d[1]) maps errors to syndromes (H_X).
D2 (C.d[2]) maps stabilizer generators to errors (H_Z^T).

Returns `true` if min( syn_block_wt / dist_block(e, im(D2)) ) >= threshold.
"""
function verify_relative_modular_expansion(C::ChainComplex, threshold::Float64; 
                                           max_module_wt::Int=2, exact::Bool=true, confidence::Float64=0.99)
    D1 = C.d[1] # H_X: C_1 -> C_0
    D2 = C.d[2] # H_Z^T: C_2 -> C_1
    
    R = base_ring(D1)
    G_grp = group(R)
    group_elements = collect(G_grp)
    lift_factor = length(G_grp)
    
    n_errors = ncols(D1) 
    
    # Lift D2 to a flat binary matrix ONCE for fast classical distance calculations
    D2_flat = _lift_matrix(D2)
    # We wrap the stabilizer generator matrix into a classical LinearCode struct 
    # so we can pass it natively into the existing distance solvers.
    StabilizerCode = LinearCode(D2_flat)
    
    current_error = [zero(R) for _ in 1:n_errors]
    
    return _relative_expansion_dfs!(
        1, 0, max_module_wt, current_error, 
        D1, StabilizerCode, R, group_elements, n_errors, threshold, lift_factor, exact, confidence
    )
end

function _relative_expansion_dfs!(start_idx::Int, current_wt::Int, max_wt::Int, 
                                  current_error::Vector{T}, D1, StabilizerCode, 
                                  R, group_elements, nc::Int, threshold::Float64, 
                                  lift_factor::Int, exact::Bool, confidence::Float64) where T
    
    if current_wt > 0
        err_matrix = matrix(R, nc, 1, current_error)
        syndrome = D1 * err_matrix
        
        # 1. Calculate Syndrome Block Weight
        syn_block_wt = 0
        for r in 1:nrows(syndrome)
            if !iszero(syndrome[r, 1])
                syn_block_wt += 1
            end
        end
        
        # 2. Short-Circuit Check: 
        # If absolute expansion is safe, relative expansion is guaranteed to be safe.
        if (syn_block_wt / current_wt) < threshold
            
            # Map the group ring vector to a flat F_2 vector
            flat_error = _lift_vector(current_error, R)
            
            # 3. Route to the Coset Solver
            # We use physical Hamming distance minimization as a proxy to find the best coset representative
            _, min_flat_err = coset_distance(StabilizerCode, flat_error; exact=exact, confidence=confidence)
            
            # Map the physical F_2 vector BACK to block weight
            rel_dist = _calculate_block_weight(min_flat_err, lift_factor)
            
            # Logical errors (rel_dist > 0, syn_wt == 0) safely fail here
            if rel_dist > 0 && (syn_block_wt / rel_dist) < threshold
                return false
            end
        end
    end
    
    if current_wt == max_wt
        return true
    end
    
    for i in start_idx:nc
        for g in group_elements
            current_error[i] = R(g)
            
            if !_relative_expansion_dfs!(i + 1, current_wt + 1, max_wt, current_error, 
                                         D1, StabilizerCode, R, group_elements, nc, threshold, 
                                         lift_factor, exact, confidence)
                return false
            end
            
            current_error[i] = zero(R)
        end
    end
    
    return true
end

"""
    coset_distance(C::AbstractLinearCode, e::Vector{Int}; 
                   exact::Bool=true, confidence::Float64=0.99, max_span::Int=15, verbose::Bool=false)

Finds the minimum weight representative of the affine coset `e + c` for `c` in `C`.
Acts as the central router, evaluating the code's geometry to dispatch to the optimal engine.
"""
function coset_distance(C::AbstractLinearCode, e::Vector{Int}; 
                        exact::Bool=true, confidence::Float64=0.99, max_span::Int=15, verbose::Bool=false)
    
    # FAST PATH: Is the error already a valid codeword? (Distance 0)
    # The parity check matrix H is the nullspace of the generator matrix.
    H = parity_check_matrix(C)
    syn = (Array(H) * e) .% Int(characteristic(C.F))
    if all(iszero, syn)
        return 0, e
    end

    if !exact
        verbose && println("Routing to Probabilistic Coset Solver (ISD)...")
        # We will modify probabilistic_minimum_distance_stern to accept target_syn
        return probabilistic_minimum_distance_stern(C; confidence=confidence, target_syn=syn, verbose=verbose)
    end
    
    # ---------------------------------------------------------
    # EXACT ROUTING: Evaluate Trellis Profile
    # ---------------------------------------------------------
    H_mat = Array(H)
    verbose && println("Evaluating exact coset router...")
    
    # Profile the trellis complexity of the parity-check matrix
    best_H, best_perm, peak_E = optimize_trellis_permutation(H_mat, 10)
    
    if peak_E <= max_span
        verbose && println("Peak E ($peak_E) <= $max_span. Routing to Exact Syndrome Trellis.")
        
        # Calculate target syndrome for the permuted parity-check matrix
        e_perm = e[best_perm]
        target_syn_perm = (best_H * e_perm) .% Int(characteristic(C.F))
        
        boundaries = optimal_sectionalization(best_H, Int(order(C.F)))
        
        # We will modify _min_weight_syndrome_sectionalized to accept target_syn
        dist, min_err_perm = _min_weight_syndrome_sectionalized(best_H, boundaries, target_syn=target_syn_perm, verbose=verbose)
        
        # Invert permutation to return physical vector to original basis
        min_err_orig = zeros(Int, C.n)
        min_err_orig[invperm(best_perm)] = min_err_perm
        
        return dist, min_err_orig
    else
        verbose && println("Peak E ($peak_E) > $max_span. Trellis intractable. Routing to Affine Brouwer-Zimmermann.")
        
        # We will modify minimum_distance to accept affine_shift_row
        return minimum_distance(C, alg=:BZ, affine_shift_row=e, verbose=verbose)
    end
end

# --- Helpers to bridge Algebra and Geometry ---

# function _lift_matrix(M)
#     # Placeholder for your existing module-to-binary lifting utility
#     # (Extracts F_2[G] blocks into classical companion blocks)
# end

# function _lift_vector(e_vec::Vector{T}, R) where T
#     G_grp = group(R)
#     N = length(G_grp)
#     flat_len = length(e_vec) * N
#     flat = zeros(Int, flat_len)
    
#     group_idx = Dict(g => i for (i, g) in enumerate(collect(G_grp)))
    
#     for (block_idx, poly) in enumerate(e_vec)
#         if !iszero(poly)
#             for (coeff, g) in terms(poly) 
#                 if coeff == 1
#                     flat[(block_idx - 1) * N + group_idx[g]] = 1
#                 end
#             end
#         end
#     end
#     return flat
# end

function _calculate_block_weight(flat_vec::Vector{Int}, lift_factor::Int)
    block_wt = 0
    num_blocks = length(flat_vec) ÷ lift_factor
    for b in 1:num_blocks
        start_idx = (b - 1) * lift_factor + 1
        end_idx = b * lift_factor
        if any(!iszero, @view flat_vec[start_idx:end_idx])
            block_wt += 1
        end
    end
    return block_wt
end
