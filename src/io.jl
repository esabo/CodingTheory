# Copyright (c) 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

function _integer_coordinate_array(
    M::CTMatrixTypes, F::CTFieldTypes
)
    expanded = _additive_expansion(M, F)
    return Int[
        Int(lift(Nemo.ZZ, expanded[r, c]))
        for r in 1:nrows(expanded), c in 1:ncols(expanded)
    ]
end

function _classical_export_matrix(
    C::AbstractLinearCode, representation::Symbol
)
    representation == :generator && return generator_matrix(C)
    representation in (:parity_check, :parity) &&
        return parity_check_matrix(C)
    throw(ArgumentError(
        "Expected representation=:generator or :parity_check."))
end

"""
$(TYPEDSIGNATURES)

Return a classical code matrix as a plain `Matrix{Int}` over the prime field.
Extension-field symbols are expanded into adjacent coordinate columns.

The default representation is the generator matrix for general linear codes
and the parity-check matrix for LDPC codes.
"""
function code_matrix_array(
    C::AbstractLinearCode; representation::Symbol=:generator
)
    M = _classical_export_matrix(C, representation)
    return _integer_coordinate_array(M, C.F)
end

function code_matrix_array(
    C::AbstractLDPCCode; representation::Symbol=:parity_check
)
    M = _classical_export_matrix(C, representation)
    return _integer_coordinate_array(M, C.F)
end

function _write_integer_csv(path::AbstractString, values::Matrix{Int})
    open(path, "w") do io
        for r in axes(values, 1)
            println(io, join(view(values, r, :), ','))
        end
    end
    return path
end

"""
$(TYPEDSIGNATURES)

Return `path` after writing `code_matrix_array(C)` as a header-free numeric CSV suitable for
`numpy.loadtxt(path, delimiter=",", dtype=int)`. For LDPC codes the default
is `representation=:parity_check`.

CSV stores only one matrix and is not a complete code serialization.
"""
function write_code_csv(
    path::AbstractString, C::AbstractLinearCode;
    representation::Symbol=:generator,
)
    return _write_integer_csv(
        path, code_matrix_array(C; representation=representation))
end

function write_code_csv(
    path::AbstractString, C::AbstractLDPCCode;
    representation::Symbol=:parity_check,
)
    return _write_integer_csv(
        path, code_matrix_array(C; representation=representation))
end

function _code_export_type(path::AbstractString, type::Symbol)
    selected = type == :nz ? :npz : type
    selected != :auto && return selected
    extension = lowercase(splitext(path)[2])
    extension == ".csv" && return :csv
    extension == ".toml" && return :toml
    extension == ".jld2" && return :jld2
    extension in (".npz", ".nz") && return :npz
    extension in (".pauli", ".stab") && return :pauli
    throw(ArgumentError(
        "Cannot infer an export type from `$path`; pass type=:csv, " *
        ":toml, :jld2, :npz, or :pauli."))
end

"""
$(TYPEDSIGNATURES)

Return `path` after exporting a classical, LDPC, stabilizer, or subsystem code through a unified
interface. The backend is selected by `Val(type)`; `type=:auto` infers it from
the file extension. `:nz` is accepted as an alias for `:npz`.

CSV and NPZ are numeric matrix exports. TOML and JLD2 are complete portable
quantum-code formats. Pauli output is a lossy binary quantum format.
"""
function save_code(
    path::AbstractString, C::AbstractCode;
    type::Symbol=:auto, kwargs...,
)
    selected = _code_export_type(path, type)
    return save_code(Val(selected), path, C; kwargs...)
end

function save_code(
    ::Val{F}, path::AbstractString, C::AbstractCode; kwargs...
) where F
    throw(ArgumentError(
        "Export type `$F` is unavailable for $(typeof(C))."))
end

function save_code(
    ::Val{:csv}, path::AbstractString, C::AbstractLinearCode;
    representation::Symbol=C isa AbstractLDPCCode ?
        :parity_check : :generator,
)
    return write_code_csv(path, C; representation=representation)
end
