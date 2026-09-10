module JLD2Ext

import CodingTheory
import CodingTheory: TriangularColorCode488, TriangularColorCode666, StabilizerCode, set_logicals!, set_minimum_distance! #PlanarSurfaceCode3D, PlanarSurfaceCode3D_X, ToricCode3D,
import JLD2
import JLD2: @load
import Oscar: GF, matrix

include("Quantum/misc_known_codes.jl")

function CodingTheory._save_quantum_code(
    ::Val{:jld2}, path::AbstractString, payload
)
    JLD2.jldsave(path; codingtheory_quantum_code=payload)
    return path
end

function CodingTheory._load_quantum_code(
    ::Val{:jld2}, path::AbstractString
)
    return JLD2.load(path, "codingtheory_quantum_code")
end

function CodingTheory.save_code(
    ::Val{:jld2}, path::AbstractString,
    C::CodingTheory.AbstractLinearCode;
    representation::Symbol=C isa CodingTheory.AbstractLDPCCode ?
        :parity_check : :generator,
)
    payload = Dict{String, Any}(
        "format" => "CodingTheory.matrix-export",
        "field_order" => Int(CodingTheory.order(C.F)),
        "characteristic" => Int(CodingTheory.characteristic(C.F)),
        "field_degree" => CodingTheory.degree(C.F),
        "length" => C.n,
        "representation" => String(representation),
        "matrix" => CodingTheory.code_matrix_array(
            C; representation=representation),
    )
    JLD2.jldsave(path; codingtheory_code_export=payload)
    return path
end

end
