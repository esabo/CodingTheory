# Copyright (c) 2021, 2022, 2023 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

# TODO: did we turn the Oscar banner off? ENV["NEMO_PRINT_BANNER"] = "false"

module CodingTheory

using AutoHashEquals
using CairoMakie, NetworkLayout, GraphMakie, GLMakie, WGLMakie, Graphs
import Graphs as Grphs
using Plots
# import Grphs: nv, incidence_matrix, vertices, inneighbors, add_edge!, add_vertex!
using Oscar
using JLD2
using Combinatorics
using .Threads
using GAP
using LinearAlgebra
using SparseArrays
using Random
using JuMP, GLPK
using DataStructures
using StatsBase

# TODO: no need to import things I switched to snake_case
import LinearAlgebra: tr, Adjoint, transpose, kron
import Oscar: dual, isprime, factor, transpose, order, polynomial, nrows, ncols, degree,
    isisomorphic, lift, quo, VectorSpace, dimension, extend, support, complement,
    is_regular, iscyclic, genus, density, isdegenerate, index, generators, copy, issubfield, ⊗,
    girth, generator_matrix, polynomial_ring, is_primitive, normal_subgroups, vector_space
import Oscar.Nemo: exponent_vectors
import Base: circshift, iseven, show, length, in, zeros, ⊆, /, *, ==, ∩, +, -, copy, isequal, ∘
import CairoMakie: save
import Combinatorics: powerset
import DataStructures: capacity

# TODO: don't want this here
# GAP.Packages.load("LINS");

#############################
         # types.jl
#############################

include("types.jl")
# classical types
export AbstractCode, AbstractNonadditiveCode, AbstractNonlinearCode, AbstractAdditiveCode,
    AbstractLinearCode, AbstractLDPCCode, AbstractMatrixProductCode, AbstractReedMullerCode,
    AbstractCyclicCode, AbstractBCHCode, AbstractReedSolomonCode, AbstractQuasiCyclicCode,
    AbstractGeneralizedReedSolomonCode, AbstractAlgebraicGeometryCode, WeightEnumerator
# quantum types
export AbstractSubsystemCode, AbstractSubsystemCodeCSS, AbstractStabilizerCode, AbstractStabilizerCodeCSS,
    AbstractGraphStateSubsystem, AbstractGraphStateSubsystemCSS, AbstractGraphStateStabilizer,
    AbstractGraphStateStabilizerCSS, AbstractHypergraphProductCode, AbstractEASubsystemCode,
    AbstractEASubsystemCodeCSS, AbstractEAStabilizerCode, AbstractEAStabilizerCodeCSS 
# misc
export LogicalTrait, GaugeTrait, HasLogicals, HasNoLogicals, HasGauges, HasNoGauges, copy

#############################
         # utils.jl
#############################

include("utils.jl")
export kronecker_product, Hamming_weight, weight, wt, Hamming_distance, distance,
    dist, tr, expand_matrix, symplectic_inner_product, are_symplectic_orthogonal,
    Hermitian_inner_product, Hermitian_conjugate_matrix, FpmattoJulia, is_triorthogonal,
    print_string_array, print_char_array, print_symplectic_array, pseudoinverse,
    quadratic_to_symplectic, symplectic_to_quadratic, _remove_empty, quadratic_residues,
    digits_to_int, is_basis, primitive_basis, #polynomial_basis, monomial_basis,
    normal_basis, dual_basis, complementary_basis, verify_dual_basis,
    verify_complementary_basis, are_equivalent_basis, is_self_dual_basis,
    is_primitive_basis, is_normal_basis, is_extension, poly_to_circ_matrix,
    is_regular, edge_vertex_incidence_matrix, edge_vertex_incidence_graph,
    is_valid_bipartition, extract_bipartition, is_Hermitian_self_orthogonal,
    row_supports, row_supports_symplectic
    # , _min_wt_row
    # , circ_shift
    # , lift
    # , _process_strings
    # , _Pauli_string_to_symplectic

#############################
        # cyclotomic.jl
#############################

include("Classical/cyclotomic.jl")
export ord, cyclotomic_coset, all_cyclotomic_cosets, complement_qcosets,
    qcoset_pairings, qcoset_pairings, qcoset_table, dual_qcosets

#############################
       # linearcode.jl
#############################

include("Classical/linearcode.jl")
export LinearCode, field, length, dimension, cardinality, rate, relative_distance, generator_matrix,
    parity_check_matrix, genus, Singleton_bound, number_correctable_errors, encode, syndrome, in, ⊆, ⊂,
    is_subcode, dual, Hermitian_dual, are_equivalent, is_self_dual, is_self_orthogonal, is_weakly_self_dual, ⊕,
    VectorSpace, set_minimum_distance!, permute_code, words, codewords, elements, is_MDS, is_even,
    is_doubly_even, is_triply_even, characteristic_polynomial, is_Hermitian_LCD, is_Hermitian_dual_containing,
    is_LCD, Hermitian_hull, hull, is_Hermitian_self_dual, is_dual_containing, minimum_distance_lower_bound,
    minimum_distance_upper_bound, set_distance_upper_bound!, standard_form_permutation, genus,
    is_overcomplete, are_permutation_equivalent, vector_space

#############################
    # newcodesfromold.jl
#############################

include("Classical/newcodesfromold.jl")
export code_complement, quo, quotient, /, direct_sum, ⊗, kron, tensor_product, direct_product, ×,
    product_code, extend, puncture, expurgate, augment, shorten, lengthen, u_u_plus_v,
    Plotkin_construction, subcode, juxtaposition, construction_X, construction_X3,
    u_plus_w_v_plus_w_u_plus_v_plus_w, entrywise_product_code, *, Schur_product_code, Hadamard_product_code,
    componentwise_product_code, expanded_code, subfield_subcode, trace_code, even_subcode,
    doubly_even_subcode, subcode_of_dimension_between_codes

#############################
        # LDPC/codes.jl
#############################

include("LDPC/codes.jl")
export variable_degree_distribution, check_degree_distribution,
    degree_distributions, column_bound, row_bound, bounds, density, is_regular,
    LDPCCode, degree_distributions_plot, variable_degree_polynomial,
    check_degree_polynomial, column_row_bounds, limited, regular_LDPC_code,
    shortest_cycle_ACE, shortest_cycles, computation_graph, ACE_distribution,
    ACE_spectrum, count_short_cycles, average_ACE_distribution,
    median_ACE_distribution, mode_ACE_distribution, girth

#############################
        # LDPCalgs.jl
#############################

include("LDPC/algorithms.jl")

#############################
     # LDPC/decoders.jl
#############################

include("LDPC/decoders.jl")
export Gallager_A, Gallager_B, sum_product, sum_product_box_plus, min_sum,
    find_MP_schedule, MPNoiseModel, decoder_simulation

#############################
     # LDPC/channels.jl
#############################

include("LDPC/channels.jl")

export BinaryErasureChannel, BEC, BinarySymmetricChannel, BSC, BAWGNChannel,
    BAWGNC, capacity

#############################
     # LDPC/analysis.jl
#############################

include("LDPC/analysis.jl")

export LDPCEnsemble, erasure_probability, crossover_probability, standard_deviation, variance,
    type, density_evolution, density_evolution!, plot_EXIT_chart, multiplicative_gap,
    multiplicative_gap_lower_bound, density_lower_bound, check_concentrated_degreedistribution

#############################
    # MatrixProductCode.jl
#############################

include("Classical/MatrixProductCode.jl")
export MatrixProductCode

#############################
      # ReedMuller.jl
#############################

include("Classical/ReedMuller.jl")
export order, RMr, RMm, ReedMullerCode, number_of_variables

#############################
      # cycliccode.jl
#############################

include("Classical/cycliccode.jl")
export defining_set, splitting_field, polynomial_ring, primitive_root, offset,
    design_distance, qcosets, qcosets_reps, generator_polynomial, parity_check_polynomial,
    idempotent, is_primitive, is_narrowsense, is_reversible, find_delta, dual_defining_set,
    CyclicCode, BCHCode, ReedSolomonCode, complement, ==, ∩, +, QuadraticResidueCode,
    zeros, BCHbound, is_degenerate, nonzeros, is_cyclic, is_antiprimitive

#############################
    # quasicycliccode.jl
#############################

include("Classical/quasicycliccode.jl")
export weight_matrix, base_matrix, protograph_matrix, QuasiCyclicCode, index,
    expansion_factor, type, polynomial_matrix, polynomial_matrix_type,
    noncirculant_generator_matrix, noncirculant_parity_check_matrix, generators,
    circulants, is_single_generator, lift

#############################
      # miscknowncodes.jl
#############################

include("Classical/miscknowncodes.jl")
export ZeroCode, IdentityCode, RepetitionCode, SingleParityCheckCode, SPCCode,
    Hexacode, HammingCode, TetraCode, SimplexCode, GolayCode, ExtendedGolayCode

#############################
 # GeneralizedReedSolomon.jl
#############################

include("Classical/GeneralizedReedSolomon.jl")
export GeneralizedReedSolomonCode, scalars, dual_scalars, evaluation_points

#############################
      # concatenation.jl
#############################

include("Classical/concatenation.jl")
export concatenate, inner_code, outer_code, expansion_basis, expansion_dual_basis, concatenation_type

#############################
      # subsystemcode.jl
#############################

include("Quantum/subsystemcode.jl")
export SubsystemCode, field, length, num_qubits, dimension, cardinality,
    rate, signs, X_signs, Z_signs, stabilizers, symplectic_stabilizers, X_stabilizers, Z_stabilizers,
    num_X_stabs, num_Z_stabs, character_vector, is_over_complete, is_CSS, relative_distance, logicals,
    logical_operators, bare_logicals, bare, logicals_matrix, gauges, gauge_operators, gauges_matrix,
    gauge_operators_matrix, dressed, dressed_operators, dressed_logicals, gauge_group, gauge_group_matrix,
    gauge_generators_matrix, gauge_group_generators_matrix, set_signs!, set_logicals!, set_minimum_distance!,
    split_stabilizers, is_logical, syndrome, X_syndrome, Z_syndrome, promote_logicals_to_gauge!, swap_X_Z_logicals!,
    swap_X_Z_gauge_operators!, all_stabilizers, elements, print_all_stabilizers, print_all_elements,
    augment, expurgate, fix_gauge, set_X_stabilizers, set_Z_stabilizers, set_stabilizers,
    set_Z_stabilizers!, set_distance_lower_bound!, permute_code!, permute_code, set_stabilizers!,
    set_X_stabilizers!, standard_form_A, standard_form_A1, standard_form_A2, standard_form_B, standard_form_C1,
    standard_form_C2, standard_form_D, standard_form_E, logicals_standard_form, promote_gauges_to_logical!,
    promote_gauges_to_logical, promote_logicals_to_gauge

#############################
      # stabilizercode.jl
#############################

include("Quantum/stabilizercode.jl")
export StabilizerCodeCSS, CSSCode, StabilizerCode

#############################
      # graphstate.jl
#############################

include("Quantum/graphstate.jl")
export ClusterState, GraphState

#############################
  # miscknownquantumcodes.jl
#############################

include("Quantum/miscknowncodes.jl")
# subsystem
export GaugedShorCode, Q9143, BaconShorCode, BravyiBaconShorCode, GeneralizedBaconShorCode,
    NappPreskill3DCode, NappPreskill4DCode, SubsystemToricCode, SubsystemSurfaceCode

# stabilizer
export FiveQubitCode, Q513, SteaneCode, Q713, _SteaneCodeTrellis, ShorCode, Q913,
    Q412, Q422, Q511, Q823, Q15RM, Q1513, Q1573, TriangularSurfaceCode,
    RotatedSurfaceCode, XZZXSurfaceCode, TriangularColorCode488, TriangularColorCode666,
    ToricCode, PlanarSurfaceCode, XYSurfaceCode, XYZ2Code, HCode, QC6, QC4, ToricCode4D

#############################
        # trellis.jl
#############################

include("trellis.jl")
export Trellis, vertices, edges, isisomorphic, isequal, loadbalancedecode,
    trellisorientedformlinear,trellisprofiles, syndrometrellis,
    trellisorientedformadditive, optimalsectionalizationQ, weightQ!,
    shiftandweightQ!, shiftanddecodeQ!, shift!, isshifted

#############################
      # weight_dist.jl
#############################

include("weight_dist.jl")
export polynomial, type, CWE_to_HWE, weight_enumerator, MacWilliams_identity,
    weight_distribution, weight_plot, support, minimum_distance, weight_plot_CSS_X,
    weight_plot_CSS_Z, weight_plot_CSS, minimum_distance_X_Z, minimum_distance_X,
    minimum_distance_Z, is_pure, Sterns_attack, Gray_code_minimum_distance, minimum_words,
    words_of_weight, QDistRndCSS

#############################
   # quantumproductcodes.jl
#############################

include("Quantum/productcodes.jl")
export HypergraphProductCode, GeneralizedShorCode, BaconCasaccinoConstruction,
    HyperBicycleCodeCSS, HyperBicycleCode, GeneralizedBicycleCode,
    GeneralizedHypergraphProductCodeMatrices, LiftedGeneralizedHypergraphProductCode,
    QuasiCyclicLiftedProductCodeMatrices, QuasiCyclicLiftedProductCode,
    BiasTailoredQuasiCyclicLiftedProductCodeMatrices, BiasTailoredQuasiCyclicLiftedProductCode

#############################
        # tilings.jl
#############################

include("tilings.jl")
export ReflectionGroup, triangle_group, r_s_group, tetrahedron_group, q_r_s_group,
    star_tetrahedron_group, cycle_tetrahedron_group, normal_subgroups, is_fixed_point_free,
    is_orientable, is_k_colorable, coset_intersection

#############################
        # Tanner.jl
#############################

include("Classical/Tanner.jl")
export Tanner_graph_plot, Tanner_graph, Tanner_code

end
