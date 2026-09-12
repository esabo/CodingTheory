module NPZExt

import CodingTheory
import NPZ

function _common_metadata(C, kind::Int)
    return Dict{String, Any}(
        "field_order" => reshape([Int(CodingTheory.order(C.F))], 1, 1),
        "characteristic" =>
            reshape([Int(CodingTheory.characteristic(C.F))], 1, 1),
        "field_degree" => reshape([CodingTheory.degree(C.F)], 1, 1),
        "length" => reshape([C.n], 1, 1),
        "code_kind" => reshape([kind], 1, 1),
    )
end

function CodingTheory.save_code(
    ::Val{:npz}, path::AbstractString,
    C::CodingTheory.AbstractLinearCode;
    representation::Symbol=C isa CodingTheory.AbstractLDPCCode ?
        :parity_check : :generator,
)
    kind = C isa CodingTheory.AbstractLDPCCode ? 1 : 0
    arrays = _common_metadata(C, kind)
    arrays["matrix"] = CodingTheory.code_matrix_array(
        C; representation=representation)
    arrays["representation"] = reshape(
        [representation in (:parity, :parity_check) ? 1 : 0], 1, 1)
    NPZ.npzwrite(path, arrays)
    return path
end

function CodingTheory.save_code(
    ::Val{:npz}, path::AbstractString,
    S::CodingTheory.AbstractSubsystemCode;
    generators::Symbol=:presentation,
)
    kind = CodingTheory.GaugeTrait(typeof(S)) ==
        CodingTheory.HasGauges() ? 3 : 2
    arrays = _common_metadata(S, kind)
    arrays["matrix"] = CodingTheory.quantum_generator_array(
        S; generators=generators)
    char_vec = CodingTheory.character_vector(S)
    isempty(char_vec) || (arrays["character_vector"] = reshape(
        Int[Int(CodingTheory.lift(CodingTheory.Nemo.ZZ, x))
            for x in char_vec], 1, :))
    NPZ.npzwrite(path, arrays)
    return path
end

end
