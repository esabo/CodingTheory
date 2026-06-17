# Copyright (c) 2021 - 2026 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

module CodingTheory

using AutoHashEquals
using Graphs
import Graphs as Grphs
using Oscar
using Combinatorics
using .Threads
using LinearAlgebra
using SparseArrays
using Random
using DataStructures
using StatsBase
using Distributions
using ProgressMeter
using DocStringExtensions
using QuadGK
using SpecialFunctions

import LinearAlgebra: tr, Adjoint, transpose, kron, diagm, dot
import Oscar: dual, factor, transpose, order, polynomial, nrows, ncols, degree,
    lift, quo, vector_space, dimension, extend, support, complement,
    is_regular, is_cyclic, genus, density, is_degenerate, index, generators, copy, is_subfield, ⊗,
    girth, generator_matrix, polynomial_ring, is_primitive, normal_subgroups, vector_space,
    tensor_product, gens, dim, is_isomorphic, field, is_irreducible, SMat, extension_field, ⊕,
    number_of_variables
import Oscar.Hecke: is_separable, ⊗, ⊕
import Oscar.Nemo: exponent_vectors
import Oscar.GAP: GapObj, Globals, Packages
import Base: circshift, iseven, show, length, in, zeros, ⊆, /, *, ==, ∩, +, -, copy, isequal, ∘, ∈,
    getproperty, setproperty!
import Combinatorics: powerset
import DataStructures: capacity
import SpecialFunctions: erfc

# tilings.jl
LINS_flag_install = Packages.install("LINS")
if LINS_flag_install
    LINS_flag = Packages.load("LINS")
    LINS_flag || @warn "Unable to load the GAP packages LINS."
else
    @warn "Unable to install the GAP packages LINS."
end

#############################
         # types.jl
#############################

const CTFieldTypes = FinField
const CTFieldElem = FinFieldElem
const CTMatrixTypes = Union{fpMatrix, FqMatrix, SparseMatrixCSC, SMat} # MatElem{<:CTFieldElem}
const CTPolyRing = PolyRing{<:CTFieldElem}
const CTPolyRingElem = PolyRingElem{<:CTFieldElem}
const CTGroupAlgebra = GroupAlgebraElem{fpFieldElem, GroupAlgebra{fpFieldElem, FinGenAbGroup, FinGenAbGroupElem}}
const CTChainComplex = Union{ComplexOfMorphisms{AbstractAlgebra.FPModule{fpFieldElem}}} # residue and group algebras later
const CTPolyMatrix = Union{AbstractAlgebra.Generic.MatSpaceElem{fpPolyRingElem}, AbstractAlgebra.Generic.MatSpaceElem{FqPolyRingElem}}

include("Classical/types.jl")
# Export Abstract Types
export AbstractCode, AbstractNonadditiveCode, AbstractNonlinearCode, AbstractAdditiveCode, 
       AbstractLinearCode, AbstractMatrixProductCode, AbstractReedMullerCode, AbstractCyclicCode, 
       AbstractBCHCode, AbstractReedSolomonCode, AbstractCyclicCode2D, AbstractQuasiCyclicCode, 
       AbstractGeneralizedReedSolomonCode, AbstractAlgebraicGeometryCode, AbstractConcatenatedCode, 
       AbstractAlternateCode, AbstractGoppaCode, AbstractGeneralizedSrivastavaCode, 
       AbstractTwistedReedSolomonCode, AbstractTannerCode

# # Export Concrete Types & Structs
# export HammingWeightEnumerator, CompleteWeightEnumerator, ExtendedQRCode, ProductCode, 
#        LinearCode, CyclicCode, BCHCode, ReedSolomonCode, QuasiCyclicCode, 
#        GeneralizedReedSolomonCode, AlternateCode, GeneralizedSrivastavaCode, ConcatenatedCode, 
#        TwistedReedSolomonCode, HammingCode, SimplexCode, MacDonaldCode, PlotkinCode, 
#        DirectSumCode, TensorProductCode, MultilevelConcatenatedCode, GabidulinCode, 
#        GoppaCode, MatrixProductCode, ReedMullerCode, TannerCode

include("LDPC/types.jl")
# Abstract Types
export AbstractLDPCCode, AbstractChannel, AbstractDiscreteChannel, 
       AbstractContinuousChannel

# Concrete Types & Aliases
export LDPCCode, BinaryErasureChannel, BEC, BinarySymmetricChannel, BSC, 
       BAWGNChannel, BAWGNC, ZChannel, RayleighFadingChannel, LDPCEnsemble, 
       METEnsemble, AbstractLDPCFamily

# include("Quantum/types.jl")
# export AbstractSubsystemCode, AbstractSubsystemCodeCSS, AbstractStabilizerCode, AbstractStabilizerCodeCSS,
#     AbstractGraphStateSubsystem, AbstractGraphStateSubsystemCSS, AbstractGraphStateStabilizer,
#     AbstractGraphStateStabilizerCSS, AbstractHypergraphProductCode, AbstractEASubsystemCode,
#     AbstractEASubsystemCodeCSS, AbstractEAStabilizerCode, AbstractEAStabilizerCodeCSS 
# # misc
# export LogicalTrait, GaugeTrait, CSSTrait, HasLogicals, HasNoLogicals, HasGauges, HasNoGauges,
#     IsCSS, IsNotCSS, copy, ChainComplex

#############################
         # utils.jl
#############################

include("utils.jl")
export kronecker_product, Hamming_weight, weight, wt, Hamming_distance, distance,
    dist, tr, expand_matrix, symplectic_inner_product, are_symplectic_orthogonal,
    Hermitian_inner_product, Hermitian_conjugate_matrix, is_triorthogonal,
    print_string_array, print_char_array, print_symplectic_array, pseudoinverse,
    quadratic_to_symplectic, symplectic_to_quadratic, _remove_empty, quadratic_residues,
    digits_to_int, is_basis, primitive_basis, #polynomial_basis, monomial_basis,
    normal_basis, dual_basis, complementary_basis, verify_dual_basis,
    verify_complementary_basis, are_equivalent_basis, is_self_dual_basis,
    is_primitive_basis, is_normal_basis, is_extension,
    is_regular, edge_vertex_incidence_matrix, edge_vertex_incidence_graph,
    is_valid_bipartition, extract_bipartition, is_Hermitian_self_orthogonal,
    row_supports, row_supports_symplectic, strongly_lower_triangular_reduction,
    residue_polynomial_to_circulant_matrix, group_algebra_element_to_circulant_matrix,
    load_alist, extended_binomial
    # load_alist, _rref_non_pivot_cols
    # , _min_wt_row
    # , circ_shift
    # , lift
    # , _process_strings
    # , _Pauli_string_to_symplectic

#############################
        # iterators.jl
#############################

# include("iterators.jl")

#############################
# Classical/concatenation.jl
#############################

include("Classical/concatenation.jl")
export concatenate, multilevel_concatenation, inner_code, outer_code, 
       expansion_basis, expansion_dual_basis, concatenation_type

#############################
  # Classical/cyclic_code.jl
#############################

include("Classical/cyclic_code.jl")
# Constructors
export CyclicCode, BCHCode, ReedSolomonCode, QuadraticResidueCode, FireCode, 
       PolyadicCodes, DuadicCodes, TriadicCodes, TetradicCodes

# Getters & Bounds
export splitting_field, polynomial_ring, primitive_root, offset, design_distance, 
       qcosets, qcosets_reps, defining_set, zeros, nonzeros, generator_polynomial, 
       parity_check_polynomial, idempotent, BCH_bound, BCH_offset, HT_bound, 
       Roos_bound, dual_defining_set, trace_representation

# Properties & Classifications
export is_cyclic, is_subcode, is_self_dual, is_narrowsense, is_reversible, 
       is_degenerate, is_primitive, is_antiprimitive, is_irreducible

# Operations & Multipliers
export complement, dual, Schur_product_code, 
       Hadamard_product_code, componentwise_product_code, apply_multiplier, 
       is_multiplier_equivalent, multiplier_group, multiplier_subgroup_Sn, 
       multiplier_subgroup_Zn, constituents, ambient_constituents, 
       MattsonSolomon_transform, inverse_MattsonSolomon_transform

#############################
  # Classical/cyclotomic.jl
#############################

include("Classical/cyclotomic.jl")
export ord, cyclotomic_coset, all_cyclotomic_cosets, complement_qcosets, 
       qcoset_pairings, qcoset_table, dual_qcosets, minimal_polynomial, 
       are_conjugates


#############################
   # Classical/designs.jl
#############################

include("Classical/designs.jl")
export is_design_holder, design_strength, minimum_weight_blocks

#############################
   # Classical/Gabidulin.jl
#############################

include("Classical/Gabidulin.jl")
export GeneralizedGabidulinCode, GabidulinCode, RandomGabidulinCode

#############################
   # Classical/Gleason.jl
#############################

include("Classical/Gleason.jl")
export Gleason_bound, is_extremal

#############################
    # Classical/Goppa.jl
#############################

include("Classical/Goppa.jl")
export GoppaCode, RandomGoppaCode, Goppa_polynomial, extension_field, 
       is_irreducible, is_seperable, nonzeros, is_cumulative

#############################
# Classical/GRS_alternate.jl
#############################

include("Classical/GRS_alternate.jl")
export GeneralizedReedSolomonCode, RandomGeneralizedReedSolomonCode, 
       GeneralizedSrivastavaCode, SrivastavaCode, GeneralizedBCHCode, 
       scalars, dual_scalars, evaluation_points, is_primitive, 
       syndromes, syndrome_polynomial, AlternateCode

#############################
# Classical/invariant_theory.jl
#############################

include("Classical/invariant_theory.jl")
export Mallows_Sloane_bound, Bachoc_Gaborit_bound, gleason_generators, 
       shadow_transform, is_valid_self_dual_enumerator, 
       extremal_weight_enumerator

#############################
 # Classical/ISD_attacks.jl
#############################

include("Classical/ISD_attacks.jl")
export Stern_attack, Prange_attack, Lee_Brickell_attack, Leon_attack, 
       Canteaut_Chabaud_attack, DOOM_Stern_attack, MMT_attack, BJMM_attack, 
       required_ISD_iterations, Gilbert_Varshamov_bound, syndrome_decode

#############################
  # Classical/linear_code.jl
#############################

include("Classical/linear_code.jl")
# Constructors & Core Structure
export LinearCode, random_linear_code, field, length, dimension, cardinality, 
       rate, generator_matrix, parity_check_matrix, standard_form_permutation,
       ⊕

# Bounds & Properties
export relative_distance, genus, minimum_distance_lower_bound, 
       minimum_distance_upper_bound, is_MDS, number_correctable_errors, 
       is_overcomplete, set_distance_lower_bound!, set_distance_upper_bound!, 
       set_minimum_distance!, change_field!, change_field, ⊆, ⊂

# Encoding & Code Space
export encode, syndrome, information_set, random_information_set, Singleton_bound, 
       vector_space, is_even, is_doubly_even, is_triply_even, words, codewords, 
       elements, in

# Duals, Hulls, & Equivalences
export dual, Euclidean_dual, Hermitian_dual, l_Galois_dual, hull, Euclidean_hull, 
       Hermitian_hull, l_Galois_hull, are_equivalent, is_self_dual, 
       is_self_orthogonal, is_weakly_self_dual, is_Euclidean_self_orthogonal, 
       is_dual_containing, is_Euclidean_dual_containing, is_Hermitian_self_dual, 
       is_Hermitian_self_orthogonal, is_Hermitian_weakly_self_dual, 
       is_Hermitian_dual_containing, is_l_Galois_self_dual, 
       is_l_Galois_self_orthogonal, is_l_Galois_weakly_self_dual, 
       is_l_Galois_dual_containing, is_LCD, is_Hermitian_LCD, is_l_Galois_LCD, 
       characteristic_polynomial, contains_self_dual_subcode

#############################
# Classical/MatrixProductCode.jl
#############################

include("Classical/MatrixProductCode.jl")
export RandomMatrixProductCode, constituent_codes, defining_matrix, MatrixProductCode

#############################
# Classical/McEliece.jl
#############################

# include("Classical/McEliece.jl")
# export McEliece_attack, generalized_Stern_attack

#############################
# Classical/min_dist_exact.jl
#############################

include("Classical/min_dist_exact.jl")
export information_sets, heuristic_info_set_selection, minimum_distance, 
       minimum_distance_zssmp

#############################
# Classical/min_dist_heuristics.jl
#############################

include("Classical/min_dist_heuristics.jl")
export heuristic_minimum_distance_ga, heuristic_minimum_distance_aco, 
       heuristic_minimum_distance_gga_order, heuristic_minimum_distance_irons, 
       heuristic_minimum_distance_nncs, heuristic_weight_distribution_hirotomo

#############################
# Classical/min_dist_probabilistic.jl
#############################

include("Classical/min_dist_probabilistic.jl")
export logbinomial, probabilistic_minimum_distance_prange, 
       probabilistic_minimum_distance_lee_brickell, probabilistic_minimum_distance_leon, 
       probabilistic_minimum_distance_stern, probabilistic_minimum_distance_stern_DOOM, 
       probabilistic_minimum_distance_mmt, probabilistic_minimum_distance_bjmm

#############################
# Classical/misc_known_codes.jl
#############################

include("Classical/misc_known_codes.jl")
export ZeroCode, IdentityCode, RepetitionCode, SingleParityCheckCode, SPCCode, 
       Hexacode, HammingCode, ExtendedHammingCode, TetraCode, SimplexCode, 
       ExtendedGolayCode, GolayCode, HadamardCode, WalshHadamardCode, WalshCode, 
       best_known_linear_code, MacDonaldCode, Lexicode

#############################
# Classical/new_codes_from_old.jl
#############################

include("Classical/new_codes_from_old.jl")
export u_u_plus_v, Plotkin_construction, u_plus_w_v_plus_w_u_plus_v_plus_w, 
       construction_A, residue_code, construction_B, construction_Y1, construction_Y, 
       construction_B2, construction_X, construction_X3, direct_sum, direct_product, 
       product_code, kron, tensor_product, entrywise_product_code, Schur_product_code, 
       Hadamard_product_code, componentwise_product_code, quo, quotient, 
       code_complement, juxtaposition, permute_code, extend, even_extension, 
       puncture, expurgate, shorten, augment, lengthen, subcode, 
       subcode_of_dimension_between_codes, expanded_code, subfield_subcode, 
       trace_code, even_subcode, doubly_even_subcode, triply_even_subcode, ⊕, ⊕, ×, /

#############################
# Classical/quasi-cyclic_code.jl
#############################

include("Classical/quasi-cyclic_code.jl")
export index, expansion_factor, type, polynomial_matrix, polynomial_matrix_type, 
       is_single_generator, weight_matrix, base_matrix, protograph_matrix, 
       noncirculant_generator_matrix, noncirculant_parity_check_matrix, generators, 
       circulants, shift_matrix, exponent_matrix, has_algebraic_4_cycle, 
       has_algebraic_cycle, component_matrices, algebraic_dimension, QuasiCyclicCode

#############################
  # Classical/ReedMuller.jl
#############################

include("Classical/ReedMuller.jl")
export ReedMullerCode, RandomPermutedReedMullerCode, RandomBooleanFunction, 
       RandomRMCoset, order, RM_r, number_of_variables, RM_m

#############################
    # Classical/Tanner.jl
#############################

include("Classical/Tanner.jl")
export Tanner_graph_plot, Tanner_graph, TannerCode, graph_eigenvalues, 
       spectral_gap, Sipser_Spielman_bound

#############################
    # Classical/trellis.jl
#############################

include("Classical/trellis.jl")
export past_future_profiles, vertex_counts, edge_counts, 
       optimize_trellis_permutation, optimal_sectionalization, Krawtchouk, 
       MacWilliams_HWE_transform

#############################
# Classical/TwistedReedSolomon.jl
#############################

include("Classical/TwistedReedSolomon.jl")
export TwistedReedSolomonCode, RandomTwistedReedSolomonCode, twist_vector, 
       hook_vector, coefficient_vector, number_of_twists


#############################
    # Classical/weight_reduction.jl
#############################

include("Classical/weight_reduction.jl")
export weight_reduction


#############################
    # Classical/words_of_weight.jl
#############################

include("Classical/words_of_weight.jl")
export words_of_minimum_weight, words_of_weight, polynomial, 
       HammingWeightEnumerator, MacWilliams_transform, weight_distribution, 
       complete_weight_distribution, weight_enumerator, complete_weight_enumerator, 
       weight_distribution_array, weight_plot

#############################
    # LDPC/algorithms.jl
#############################

include("LDPC/algorithms.jl")
export PEGWorkspace, generate_peg, generate_qc_peg, generate_protograph_peg, 
       generate_mackay_neal, generate_eg2_ldpc, generate_spatially_coupled, 
       generate_eg_ldpc, generate_pg_ldpc

#############################
     # LDPC/analysis.jl
#############################

include("LDPC/analysis.jl")
export LDPCEnsemble, density_evolution, multiplicative_gap, 
       multiplicative_gap_lower_bound, density_lower_bound, 
       check_concentrated_degree_distribution, optimal_lambda, optimal_rho, 
       optimal_lambda_and_rho, optimal_threshold, EXIT_chart_data, 
       EXIT_chart_plot, protograph_threshold, PEXIT_chart_data, 
       density_evolution_MET_GA

#############################
     # LDPC/channels.jl
#############################

include("LDPC/channels.jl")
export erasure_probability, crossover_probability, standard_deviation, 
       variance, capacity, transmit, llr

#############################
        # LDPC/codes.jl
#############################

include("LDPC/codes.jl")
export LDPCCode, regular_LDPC_code, variable_degree_distribution, 
       check_degree_distribution, degree_distributions, column_bound, 
       row_bound, column_row_bounds, limited, density, is_regular, 
       variable_degree_polynomial, check_degree_polynomial, dimension, 
       design_dimension, design_rate, rate, degree_distributions_plot

#############################
     # LDPC/cycles.jl
#############################

include("LDPC/cycles.jl")
export node_adjacencies, girth, local_girth, computation_graph, remove_cycles, 
       enumerate_simple_cycles, simple_cycle_length_distribution, 
       simple_cycle_length_distribution_plot, average_simple_cycle_length, 
       median_simple_cycle_length, mode_simple_cycle_length, count_simple_cycles, 
       simple_cycle_distribution_by_variable_node, 
       simple_cycle_distribution_by_variable_node_plot, enumerate_short_cycles, 
       short_cycle_length_distribution, short_cycle_length_distribution_plot, 
       average_short_cycle_length, median_short_cycle_length, mode_short_cycle_length, 
       count_short_cycles, short_cycle_distribution_by_variable_node, 
       short_cycle_distribution_by_variable_node_plot, ACE_spectrum, 
       ACE_spectrum_plot, ACE_distribution, average_ACE_distribution, 
       median_ACE_distribution, mode_ACE_distribution

#############################
     # LDPC/decoder_post.jl
#############################

include("LDPC/decoder_post.jl")
export OSDWorkspace, init_osd_workspace, osd_decode!, 
       GRANDWorkspace, init_grand_workspace, grand_decode!, 
       WBFWorkspace, init_wbf_workspace, wbf_decode!

#############################
        # LDPC/GBP.jl
#############################

include("LDPC/GBP.jl")
export Region, RegionGraph, id, parents, ancestors, subregions, overcounting_number, 
       regions, base_regions, leaves, outer_regions, basic_clusters, 
       canonical_region_graph, region_graph_from_base_nodes, is_valid_region_graph, 
       remove_zero_overcounting_numbers, remove_generational_skips, 
       triangulate_base_regions, message_passing_order, GBPWorkspace, 
       init_gbp_workspace, gbp_decode!, init_region_beliefs!, extract_hard_decisions

#############################
    # LDPC/LP_decoders.jl
#############################

include("LDPC/LP_decoders.jl")
export LP_decoder_LDPC

#############################
    # LDPC/MP_decoders.jl
#############################

include("LDPC/MP_decoders.jl")
export DecoderWorkspace, HardDecisionWorkspace, SoftDecisionWorkspace, 
       init_hard_workspace, load_hard_channel!, init_soft_workspace, 
       load_soft_channel!, decode!, boxplus_exact, boxplus_minsum, 
       boxplus_minsum_correction, layered_schedule, balance_of_layered_schedule

#############################
    # LDPC/simulations.jl
#############################

# include("LDPC/simulations.jl")


   











# #############################
#  # Quantum/subsystem_code.jl
# #############################

# include("Quantum/subsystem_code.jl")
# export SubsystemCode, field, length, num_qubits, dimension, cardinality,
#     rate, signs, X_signs, Z_signs, stabilizers, symplectic_stabilizers, X_stabilizers, Z_stabilizers,
#     num_X_stabs, num_Z_stabs, character_vector, is_over_complete, is_CSS, relative_distance, logicals,
#     logical_operators, bare_logicals, bare, logicals_matrix, gauges, gauge_operators, gauges_matrix,
#     gauge_operators_matrix, dressed, dressed_operators, dressed_logicals, gauge_group, gauge_group_matrix,
#     gauge_generators_matrix, gauge_group_generators_matrix, set_signs!, set_logicals!, set_minimum_distance!,
#     split_stabilizers, is_logical, syndrome, X_syndrome, Z_syndrome, promote_logicals_to_gauge!, swap_X_Z_logicals!,
#     swap_X_Z_gauge_operators!, all_stabilizers, elements, print_all_stabilizers, print_all_elements,
#     augment, expurgate, fix_gauge, set_X_stabilizers, set_Z_stabilizers, set_stabilizers,
#     set_Z_stabilizers!, set_distance_lower_bound!, permute_code!, permute_code, set_stabilizers!,
#     set_X_stabilizers!, standard_form_A, standard_form_A1, standard_form_A2, standard_form_B, standard_form_C1,
#     standard_form_C2, standard_form_D, standard_form_E, logicals_standard_form, promote_gauges_to_logical!,
#     promote_gauges_to_logical, promote_logicals_to_gauge, bare_minimum_distance_lower_bound,
#     bare_minimum_distance_upper_bound, dressed_minimum_distance_lower_bound,
#     dressed_minimum_distance_upper_bound, bare_X_minimum_distance_lower_bound,
#     bare_X_minimum_distance_upper_bound, dressed_X_minimum_distance_lower_bound,
#     bare_Z_minimum_distance_lower_bound, bare_Z_minimum_distance_upper_bound,
#     dressed_Z_minimum_distance_lower_bound, dressed_Z_minimum_distance_upper_bound,
#     X_minimum_distance, Z_minimum_distance, XZ_minimum_distance, set_bare_minimum_distance!,
#     set_bare_X_minimum_distance!, set_bare_Z_minimum_distance!, set_dressed_minimum_distance!,
#     set_dressed_X_minimum_distance!, set_dressed_Z_minimum_distance!

# #############################
#  # Quantum/stabilizer_code.jl
# #############################

# include("Quantum/stabilizer_code.jl")
# export StabilizerCodeCSS, CSSCode, StabilizerCode, random_CSS_code, is_CSS_T_code,
#     minimum_distance_lower_bound, minimum_distance_upper_bound, X_minimum_distance_lower_bound,
#     X_minimum_distance_upper_bound, Z_minimum_distance_lower_bound, Z_minimum_distance_upper_bound,
#     set_X_minimum_distance!, set_Z_minimum_distance!

# #############################
#    # Quantum/graphstate.jl
# #############################

# include("Quantum/graph_state.jl")
# export ClusterState, GraphState

# #############################
# # Quantum/misc_known_codes.jl
# #############################

# include("Quantum/misc_known_codes.jl")
# # subsystem
# export GaugedShorCode, Q9143, BaconShorCode, BravyiBaconShorCode, GeneralizedBaconShorCode,
#     NappPreskill3DCode, NappPreskill4DCode, SubsystemToricCode, SubsystemSurfaceCode

# # stabilizer
# export FiveQubitCode, Q513, SteaneCode, Q713, _SteaneCodeTrellis, ShorCode, Q913,
#     Q412, Q422, Q511, Q823, Q15RM, Q1513, Q1573, TriangularSurfaceCode,
#     RotatedSurfaceCode, XZZXSurfaceCode, TriangularColorCode488, TriangularColorCode666,
#     ToricCode, PlanarSurfaceCode, XYSurfaceCode, XYZ2Code, HCode, QC6, QC4, ToricCode4D,
#     Q832, SmallestInterestingColorCode, GrossCode #, PlanarSurfaceCode3D, ToricCode3D

# #############################
#    # Quantum/weight_dist.jl
# #############################

# # include("Quantum/weight_dist.jl")
# # # export weight_plot_CSS_X, weight_plot_CSS_Z, weight_plot_CSS, minimum_distance_X_Z,
# # #     minimum_distance_X, minimum_distance_Z, is_pure, QDistRndCSS
# # export minimum_distance_upper_bound!, random_information_set_minimum_distance_bound!,
# #     QDistRnd!

# #############################
# #  Quantum/product_codes.jl
# #############################

# include("Quantum/product_codes.jl")
# export HypergraphProductCode, GeneralizedShorCode, BaconCasaccinoConstruction,
#     HyperBicycleCodeCSS, HyperBicycleCode, GeneralizedBicycleCode,
#     generalized_hypergraph_product_matrices, GHGP_matrices, lifted_product_matrices,
#     GeneralizedHypergraphProductCode, LiftedProductCode, bias_tailored_lifted_product_matrices,
#     BiasTailoredLiftedProductCode, SPCDFoldProductCode, SingleParityCheckDFoldProductCode,
#     Quintavalle_basis, asymmetric_product, symmetric_product, random_homological_product_code,
#     homological_product, ⊠, BivariateBicycleCode, CoprimeBivariateBicycleCode

# #############################
# #   Quantum/simulation.jl
# #############################

# include("Quantum/simulation.jl")
# export CSS_decoder_test, CSS_decoder_with_Bayes

# #############################
# #  Quantum/decoders/OTF.jl
# #############################

# include("Quantum/decoders/OTF.jl")
# export ordered_Tanner_forest

# #############################
#         # tilings.jl
# #############################

# include("tilings.jl")
# export ReflectionGroup, triangle_group, r_s_group, tetrahedron_group, q_r_s_group,
#     star_tetrahedron_group, cycle_tetrahedron_group, normal_subgroups, is_fixed_point_free,
#     is_orientable, is_k_colorable, coset_intersection



# #############################
#        # chaincomplex.jl
# #############################

# # put into Oscar v14
# # struct ChainComplex{T <: CTMatrixTypes}
# #     F::CTFieldTypes
# #     length::UInt8
# #     boundaries::Vector{T}
# # end

# # include("chaincomplex.jl")
# # export boundaries, cochain, distance_balancing

# #############################
# # Quantum/weight_reduction.jl
# #############################

# include("Quantum/weight_reduction.jl")
# export copying, gauging, thickening_and_choose_heights, coning, quantum_weight_reduction,
#     copying_as_coning, gauging_as_coning


# #############################
# # Quantum/homological_measurements.jl
# #############################

# include("Quantum/homological_measurements.jl")
# export homological_measurement, Cheeger_constant

end
