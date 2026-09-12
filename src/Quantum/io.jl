# Copyright (c) 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

const _QUANTUM_IO_VERSION = 1
const _QUANTUM_CACHE_KEYS = (
    :d, :dx, :dz, :d_bare, :d_dressed, :dx_bare, :dz_bare,
    :dx_dressed, :dz_dressed, :l_bound, :u_bound, :l_bound_dx,
    :u_bound_dx, :l_bound_dz, :u_bound_dz, :l_bound_bare,
    :u_bound_bare, :l_bound_dressed, :u_bound_dressed,
)

function _portable_matrix_data(M::CTMatrixTypes, F::CTFieldTypes)
    expanded = _additive_expansion(M, F)
    data = Int[
        Int(lift(Nemo.ZZ, expanded[r, c]))
        for r in 1:nrows(expanded) for c in 1:ncols(expanded)
    ]
    return Dict{String, Any}(
        "rows" => nrows(M),
        "cols" => ncols(M),
        "coordinate_degree" => degree(F),
        "data" => data,
    )
end

function _matrix_from_portable_data(
    payload::AbstractDict, F::CTFieldTypes
)
    nr = Int(payload["rows"])
    nc = Int(payload["cols"])
    m = Int(payload["coordinate_degree"])
    m == degree(F) ||
        throw(ArgumentError("Stored coordinate degree does not match the field."))
    values = payload["data"]
    values isa AbstractVector ||
        throw(ArgumentError("Stored matrix coordinates must be a flat vector."))
    length(values) == nr * nc * m ||
        throw(ArgumentError("Stored matrix has the wrong number of coordinates."))

    prime_field = _prime_subfield(F)
    p = Int(characteristic(F))
    all(x -> x isa Integer && 0 <= x < p, values) ||
        throw(ArgumentError(
            "Stored matrix coordinates must be integers in 0:$(p - 1)."))
    expanded = matrix(prime_field, nr, nc * m,
        [prime_field(Int(x)) for x in values])
    basis = m == 1 ? [one(F)] : first(primitive_basis(F, prime_field))
    result = zero_matrix(F, nr, nc)
    for r in 1:nr, c in 1:nc, j in 1:m
        coefficient = expanded[r, (c - 1) * m + j]
        iszero(coefficient) && continue
        result[r, c] += _lift_prime_element(coefficient, F) * basis[j]
    end
    return result
end

function _portable_character_vector(S::AbstractSubsystemCode)
    return Int[Int(lift(Nemo.ZZ, x)) for x in character_vector(S)]
end

function _portable_cache(S::AbstractSubsystemCode)
    result = Dict{String, Any}()
    for key in _QUANTUM_CACHE_KEYS
        value = get(S.cache, key, missing)
        value isa Integer && (result[String(key)] = Int(value))
    end
    return result
end

function _portable_weight_enumerators(S::AbstractSubsystemCode)
    haskey(S.cache, :SL_weight_enum) ||
        return Dict{String, Any}()
    enumerator = S.cache[:SL_weight_enum]
    enumerator isa ShorLaflammeWeightEnumerator ||
        return Dict{String, Any}()
    encode(hwe) = [string(get(hwe.counts, w, BigInt(0)))
                   for w in 0:hwe.n]
    return Dict{String, Any}(
        "A" => encode(enumerator.A),
        "B" => encode(enumerator.B),
    )
end

function _quantum_data_fingerprint(payload::AbstractDict)
    parts = String[]
    for key in ("format", "version", "kind", "css", "field_order",
                "characteristic", "field_degree", "length", "sparse")
        push!(parts, "$key=$(payload[key])")
    end
    matrix_data = payload["generators"]
    for key in ("rows", "cols", "coordinate_degree")
        push!(parts, "generators.$key=$(matrix_data[key])")
    end
    push!(parts, "generators.data=" *
        join(string.(matrix_data["data"]), ","))
    push!(parts, "character_vector=" *
        join(string.(payload["character_vector"]), ","))
    for (key, value) in sort!(collect(payload["cache"]); by=first)
        push!(parts, "cache.$key=$value")
    end
    enumerators = payload["weight_enumerators"]
    for key in sort!(collect(keys(enumerators)))
        push!(parts, "weight_enumerators.$key=" *
            join(string.(enumerators[key]), ","))
    end
    return bytes2hex(SHA.sha256(join(parts, '\n')))
end

"""
$(TYPEDSIGNATURES)

Return a language-neutral dictionary describing `S`. Matrices are flattened in
row-major order after expanding each `GF(q)` entry over its prime field. This
is the schema used by `save_quantum_code` and is directly consumable from
Python, Julia, or other TOML/JLD2 readers.
"""
function quantum_code_data(S::AbstractSubsystemCode)
    has_gauges = GaugeTrait(typeof(S)) == HasGauges()
    generators = has_gauges ? gauge_group(S) : stabilizers(S)
    payload = Dict{String, Any}(
        "format" => "CodingTheory.quantum",
        "version" => _QUANTUM_IO_VERSION,
        "kind" => has_gauges ? "subsystem" : "stabilizer",
        "css" => is_CSS(S),
        "field_order" => Int(order(S.F)),
        "characteristic" => Int(characteristic(S.F)),
        "field_degree" => degree(S.F),
        "length" => S.n,
        "sparse" => _is_sparse_code_matrix(generators),
        "generators" => _portable_matrix_data(generators, S.F),
        "character_vector" => _portable_character_vector(S),
        "cache" => _portable_cache(S),
        "weight_enumerators" => _portable_weight_enumerators(S),
    )
    payload["integrity_sha256"] = _quantum_data_fingerprint(payload)
    return payload
end

function _field_from_quantum_data(payload::AbstractDict)
    q = Int(payload["field_order"])
    p = Int(payload["characteristic"])
    m = Int(payload["field_degree"])
    q == p^m ||
        throw(ArgumentError("Stored field metadata is inconsistent."))
    return m == 1 ? Oscar.Nemo.Native.GF(p) : GF(q)
end

"""
$(TYPEDSIGNATURES)

Return a stabilizer or subsystem code reconstructed from a portable
`quantum_code_data` dictionary. Validate the integrity fingerprint before
restoring generators or certified cache entries.
"""
function quantum_code_from_data(
    payload::AbstractDict; restore_cache::Bool=true
)
    get(payload, "format", nothing) == "CodingTheory.quantum" ||
        throw(ArgumentError("Not a CodingTheory quantum-code payload."))
    Int(get(payload, "version", 0)) == _QUANTUM_IO_VERSION ||
        throw(ArgumentError("Unsupported quantum-code format version."))
    stored_fingerprint = get(payload, "integrity_sha256", nothing)
    stored_fingerprint === nothing &&
        throw(ArgumentError(
            "Quantum-code payload is missing its integrity fingerprint."))
    stored_fingerprint == _quantum_data_fingerprint(payload) ||
        throw(ArgumentError(
            "Quantum-code payload failed its integrity check; " *
            "generators or certified metadata were modified."))

    F = _field_from_quantum_data(payload)
    generators = _matrix_from_portable_data(payload["generators"], F)
    Bool(get(payload, "sparse", false)) &&
        (generators = _sparse_code_matrix(generators))

    p = Int(characteristic(F))
    raw_char_data = get(payload, "character_vector", Int[])
    raw_char_data isa AbstractVector ||
        throw(ArgumentError("The character vector must be a flat vector."))
    char_data = Int.(raw_char_data)
    phase_modulus = p == 2 ? 4 : p
    all(x -> 0 <= x < phase_modulus, char_data) ||
        throw(ArgumentError(
            "Character-vector entries must lie in 0:$(phase_modulus - 1)."))
    char_vec = if isempty(char_data)
        missing
    else
        R, _ = residue_ring(Nemo.ZZ, p == 2 ? 4 : p)
        [R(x) for x in char_data]
    end

    kind = String(payload["kind"])
    S = if kind == "stabilizer"
        StabilizerCode(generators; char_vec=char_vec, logs_alg=:sys_eqs)
    elseif kind == "subsystem"
        SubsystemCode(generators; char_vec=char_vec)
    else
        throw(ArgumentError("Unknown stored quantum-code kind: $kind"))
    end

    restore_cache || return S
    cache = get(payload, "cache", Dict{String, Any}())
    for (key, value) in cache
        value isa Integer && (S.cache[Symbol(key)] = Int(value))
    end
    enumerators =
        get(payload, "weight_enumerators", Dict{String, Any}())
    if haskey(enumerators, "A") && haskey(enumerators, "B")
        decode(values) = Dict(
            w - 1 => parse(BigInt, String(value))
            for (w, value) in enumerate(values) if parse(BigInt, String(value)) != 0
        )
        A = HammingWeightEnumerator(S.n, decode(enumerators["A"]))
        B = HammingWeightEnumerator(S.n, decode(enumerators["B"]))
        S.cache[:weight_enum_A] = A
        S.cache[:weight_enum_B] = B
        S.cache[:weight_dist_A] = A.counts
        S.cache[:weight_dist_B] = B.counts
        S.cache[:SL_weight_enum] = _validate_SL_enumerator(
            ShorLaflammeWeightEnumerator(S.n, A, B),
            cardinality(S), Int(order(S.F)))
    end
    return S
end

function _quantum_io_format(path::AbstractString, format::Symbol)
    format != :auto && return format
    extension = lowercase(splitext(path)[2])
    extension == ".toml" && return :toml
    extension == ".jld2" && return :jld2
    extension in (".pauli", ".stab") && return :pauli
    throw(ArgumentError(
        "Cannot infer format from `$path`; use format=:toml, :jld2, or :pauli."))
end

function _save_quantum_code(::Val{F}, path, payload) where F
    throw(ArgumentError("Quantum-code format `$F` is unavailable."))
end

function _load_quantum_code(::Val{F}, path) where F
    throw(ArgumentError("Quantum-code format `$F` is unavailable."))
end

"""
$(TYPEDSIGNATURES)

Return `path` after persisting a stabilizer or subsystem code. TOML is the portable, dependency-free
format. JLD2 is available when JLD2 is loaded. Certified cache data is bound
to the generator payload by a SHA-256 integrity fingerprint. Files ending in
`.pauli` or `.stab` omit phase/cache metadata and therefore require
`allow_lossy=true`.
"""
function save_quantum_code(
    path::AbstractString, S::AbstractSubsystemCode;
    format::Symbol=:auto, allow_lossy::Bool=false,
)
    selected = _quantum_io_format(path, format)
    if selected == :toml
        open(path, "w") do io
            TOML.print(io, quantum_code_data(S); sorted=true)
        end
    elseif selected == :pauli
        allow_lossy ||
            throw(ArgumentError(
                "Pauli files omit cache and phase metadata; " *
                "pass allow_lossy=true or use TOML/JLD2."))
        write_pauli_strings(path, S)
    else
        _save_quantum_code(Val(selected), path, quantum_code_data(S))
    end
    return path
end

"""
$(TYPEDSIGNATURES)

Return a stabilizer or subsystem code loaded from TOML, JLD2, or a phase-free
Pauli-string file. Set `restore_cache=false` to ignore stored cache entries.
"""
function load_quantum_code(
    path::AbstractString;
    format::Symbol=:auto, restore_cache::Bool=true,
)
    selected = _quantum_io_format(path, format)
    selected == :toml &&
        return quantum_code_from_data(
            TOML.parsefile(path); restore_cache=restore_cache)
    selected == :pauli && return read_pauli_strings(path)
    return quantum_code_from_data(
        _load_quantum_code(Val(selected), path);
        restore_cache=restore_cache)
end

function save_code(
    ::Val{:csv}, path::AbstractString, S::AbstractSubsystemCode;
    generators::Symbol=:presentation,
)
    return write_quantum_csv(path, S; generators=generators)
end

function save_code(
    ::Val{:toml}, path::AbstractString, S::AbstractSubsystemCode;
    kwargs...,
)
    return save_quantum_code(path, S; format=:toml, kwargs...)
end

function save_code(
    ::Val{:jld2}, path::AbstractString, S::AbstractSubsystemCode;
    kwargs...,
)
    return save_quantum_code(path, S; format=:jld2, kwargs...)
end

function save_code(
    ::Val{:pauli}, path::AbstractString, S::AbstractSubsystemCode;
    kwargs...,
)
    return save_quantum_code(path, S; format=:pauli, kwargs...)
end

function _pauli_generator_matrix(
    S::AbstractSubsystemCode, generators::Symbol
)
    if generators == :stabilizers
        return stabilizers(S)
    elseif generators == :gauges
        return gauges_matrix(S)
    elseif generators == :gauge_group
        return gauge_group(S)
    elseif generators == :logicals
        return logicals_matrix(S)
    end
    throw(ArgumentError(
        "Expected generators=:stabilizers, :gauges, :gauge_group, or :logicals."))
end

"""
$(TYPEDSIGNATURES)

Return the selected quantum generators as a plain `Matrix{Int}` over the
prime field. Extension-field entries are expanded into adjacent coordinate
columns. The result can be passed directly to PythonCall/NumPy or written with
`write_quantum_csv`.

`generators=:presentation` selects the gauge group for subsystem codes and
the stabilizer group otherwise.
"""
function quantum_generator_array(
    S::AbstractSubsystemCode; generators::Symbol=:presentation
)
    selected = if generators == :presentation
        GaugeTrait(typeof(S)) == HasGauges() ? :gauge_group : :stabilizers
    else
        generators
    end
    M = _pauli_generator_matrix(S, selected)
    return _integer_coordinate_array(M, S.F)
end

"""
$(TYPEDSIGNATURES)

Return `path` after writing a header-free numeric CSV containing `quantum_generator_array(S)`.
It is directly readable with `numpy.loadtxt(path, delimiter=",",
dtype=int)`.

CSV stores only a selected generator matrix. It does not preserve field,
phase, code-kind, or cached metadata; use TOML for a complete portable
round trip.
"""
function write_quantum_csv(
    path::AbstractString, S::AbstractSubsystemCode;
    generators::Symbol=:presentation,
)
    values = quantum_generator_array(S; generators=generators)
    return _write_integer_csv(path, values)
end

"""
$(TYPEDSIGNATURES)

Return binary symplectic generators as strings over `I`, `X`, `Y`, and `Z`.
Phase information is not included; use `save_quantum_code` for a lossless
round trip.
"""
function pauli_strings(
    S::AbstractSubsystemCode; generators::Symbol=:stabilizers
)
    Int(order(S.F)) == 2 && degree(S.F) == 1 ||
        throw(ArgumentError("Pauli strings are only defined for binary codes."))
    M = _dense_code_matrix(_pauli_generator_matrix(S, generators), S.F)
    strings = String[]
    for r in 1:nrows(M)
        chars = Vector{Char}(undef, S.n)
        for q in 1:S.n
            x = !iszero(M[r, q])
            z = !iszero(M[r, S.n + q])
            chars[q] = x ? (z ? 'Y' : 'X') : (z ? 'Z' : 'I')
        end
        push!(strings, String(chars))
    end
    return strings
end

"""
$(TYPEDSIGNATURES)

Return `path` after writing one phase-free binary Pauli generator per line.
Throw an error when the code has a nonempty character vector because this
format cannot preserve explicit phases.
"""
function write_pauli_strings(
    path::AbstractString, S::AbstractSubsystemCode;
    generators::Symbol=:stabilizers
)
    isempty(character_vector(S)) ||
        throw(ArgumentError(
            "Pauli-string output cannot preserve a nonempty character vector; " *
            "use TOML or JLD2."))
    open(path, "w") do io
        for string in pauli_strings(S; generators=generators)
            println(io, string)
        end
    end
    return path
end

"""
$(TYPEDSIGNATURES)

Return a stabilizer code, or a subsystem code when `subsystem=true`, parsed
from a file containing one phase-free binary Pauli generator per line.
"""
function read_pauli_strings(
    path::AbstractString; subsystem::Bool=false
)
    strings = String.(filter(!isempty, strip.(readlines(path))))
    isempty(strings) &&
        throw(ArgumentError("The Pauli-string file is empty."))
    any(s -> startswith(s, '+') || startswith(s, '-'), strings) &&
        throw(ArgumentError(
            "Signed Pauli strings are unsupported because their phases " *
            "cannot be represented losslessly; use TOML or JLD2."))
    length(unique(length.(strings))) == 1 ||
        throw(ArgumentError("All Pauli strings must have the same length."))
    all(s -> all(c -> c in ('I', 'X', 'Y', 'Z'), s), strings) ||
        throw(ArgumentError(
            "Pauli strings may contain only I, X, Y, and Z."))
    return subsystem ? SubsystemCode(strings) : StabilizerCode(strings)
end
