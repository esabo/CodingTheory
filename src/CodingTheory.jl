# Copyright (c) 2021, 2022, 2023 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

module CodingTheory

# TODO: make namespaces minimal
using AutoHashEquals
using Oscar
using CairoMakie, Graphs
using JLD2
using Combinatorics
using .Threads
using GAP
using LinearAlgebra # need to remove? want only mul!
using SparseArrays
using Random
using JuMP, GLPK

import LinearAlgebra: tr, Adjoint, transpose
import Oscar: dual, isprime, factor, transpose, order, polynomial, nrows, ncols, degree,
    isisomorphic, lift, quo, VectorSpace, dimension, extend, support, complement, isprimitive,
    isregular, iscyclic, genus, density, isdegenerate, index, generators, copy, issubfield
import Base: circshift, iseven, show, length, in, zeros, ⊆, /, *, ==, ∩, +, -, copy, isequal, ∘
import CairoMakie: save
import Combinatorics: powerset
import Graphs: nv, incidence_matrix

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
export kroneckerproduct, Hammingweight, weight, wt, Hammingdistance, distance,
    dist, tr, expandmatrix, symplecticinnerproduct, aresymplecticorthogonal,
    Hermitianinnerproduct, Hermitianconjugatematrix, FpmattoJulia, istriorthogonal,
    printstringarray, printchararray, printsymplecticarray, pseudoinverse,
    quadratictosymplectic, symplectictoquadratic, _removeempty, quadraticresidues,
    digitstoint, isbasis, primitivebasis, #polynomialbasis, monomialbasis,
    normalbasis, dualbasis, complementarybasis, verifydualbasis,
    verifycomplementarybasis, areequivalentbasis, isselfdualbasis,
    isprimitivebasis, isnormalbasis, isextension, polytocircmatrix,
    isregular, edgevertexincidencematrix, edgevertexincidencegraph,
    isvalidbipartition, extractbipartition, isHermitianselforthogonal,
    rowsupports, rowsupportssymplectic
    # , _minwtrow
    # , circshift
    # , lift
    # , _processstrings
    # , _Paulistringtosymplectic

#############################
        # cyclotomic.jl
#############################

include("cyclotomic.jl")
export ord, cyclotomiccoset, allcyclotomiccosets, complementqcosets,
    qcosetpairings, qcosetpairings, qcosettable, dualqcosets

#############################
       # linearcode.jl
#############################

include("linearcode.jl")
export LinearCode, field, length, dimension, cardinality, rate, relativedistance, generatormatrix,
    paritycheckmatrix, genus, Singletonbound, numbercorrectableerrors, encode, syndrome, in, ⊆, ⊂,
    issubcode, dual, Hermitiandual, areequivalent, isselfdual, isselforthogonal, isweaklyselfdual, ⊕,
    VectorSpace, setminimumdistance!, permutecode, words, codewords, elements, isMDS, iseven,
    isdoublyeven, istriplyeven, characteristicpolynomial, isHermitianLCD, isHermitiandualcontaining,
    isLCD, Hermitianhull, hull, isHermitianselfdual, isdualcontaining, minimumdistancelowerbound,
    minimumdistanceupperbound, setdistanceupperbound!, standardformpermutation, genus,
    isovercomplete, arepermutationequivalent

#############################
    # newcodesfromold.jl
#############################

include("newcodesfromold.jl")
export codecomplement, quo, quotient, /, directsum, ⊗, kron, tensorproduct, directproduct, ×,
    productcode, extend, puncture, expurgate, augment, shorten, lengthen, uuplusv,
    Plotkinconstruction, subcode, juxtaposition, constructionX, constructionX3,
    upluswvpluswuplusvplusw, entrywiseproductcode, *, Schurproductcode, Hadamardproductcode,
    componentwiseproductcode, expandedcode, subfieldsubcode, tracecode, evensubcode,
    doublyevensubcode, subcodeofdimensionbetweencodes

#############################
         # LDPC.jl
#############################

include("LDPC/LDPC.jl")
export variabledegreedistribution, checkdegreedistribution,
    degreedistributions, columnbound, rowbound, bounds, density, isregular,
    LDPCCode, degreedistributionssplot, variabledegreepolynomial,
    checkdegreepolynomial, columnrowbounds, limited, regularLDPCCode

#############################
        # LDPCalgs.jl
#############################

include("LDPC/algorithms.jl")

#############################
     # LDPC/decoders.jl
#############################

include("LDPC/decoders.jl")
export GallagerA, GallagerB, sumproduct, sumproductboxplus, minsum,
    findMPschedule, MPNoiseModel, decodersimulation

#############################
     # LDPC/analysis.jl
#############################

include("LDPC/analysis.jl")

export BinaryErasureChannel, BEC, BinarySymmetricChannel, BSC, BAWGNChannel,
    BAWGNC, LDPCEnsemble

export erasureprobability, crossoverprobability, standarddeviation, variance, capacity,
    type, densityevolution!, plotEXITchart, multiplicativegap, multiplicativegaplowerbound,
    densitylowerbound, checkconcentrateddegreedistribution, optimallambda, optimalrho,
    optimallambdaandrho, optimalthreshold

#############################
    # MatrixProductCode.jl
#############################

include("MatrixProductCode.jl")
export MatrixProductCode

#############################
      # ReedMuller.jl
#############################

include("ReedMuller.jl")
export order, RMr, RMm, ReedMullerCode, numberofvariables

#############################
      # cycliccode.jl
#############################

include("cycliccode.jl")
export definingset, splittingfield, polynomialring, primitiveroot, offset,
    designdistance, qcosets, qcosetsreps, generatorpolynomial, paritycheckpolynomial,
    idempotent, isprimitive, isnarrowsense, isreversible, finddelta, dualdefiningset,
    CyclicCode, BCHCode, ReedSolomonCode, complement, ==, ∩, +, QuadraticResidueCode,
    zeros, BCHbound, mindistlowerbound, isdegenerate, nonzeros, iscyclic, isantiprimitive,
    isdegenerate, iscyclic

#############################
    # quasicycliccode.jl
#############################

include("quasicycliccode.jl")
export weightmatrix, basematrix, protographmatrix, QuasiCyclicCode, index,
    expansionfactor, type, polynomialmatrix, polynomialmatrixtype,
    noncirculantgeneratormatrix, noncirculantparitycheckmatrix, generators,
    circulants, issinglegenerator, lift

#############################
      # miscknowncodes.jl
#############################

include("miscknowncodes.jl")
export ZeroCode, IdentityCode, RepetitionCode, SingleParityCheckCode, SPCCode,
    Hexacode, HammingCode, TetraCode, SimplexCode, GolayCode, ExtendedGolayCode

#############################
 # GeneralizedReedSolomon.jl
#############################

include("GeneralizedReedSolomon.jl")
export GeneralizedReedSolomonCode, scalars, dualscalars, evaluationpoints

#############################
      # concatenation.jl
#############################

include("concatenation.jl")
export concatenate, innercode, outercode, expansionbasis, expansiondualbasis, concatenationtype

#############################
      # subsystemcode.jl
#############################

include("subsystemcode.jl")
export SubsystemCode, field, quadraticfield, length, numqubits, dimension, cardinality,
    rate, signs, Xsigns, Zsigns, stabilizers, symplecticstabilizers, Xstabilizers, Zstabilizers,
    numXstabs, numZstabs, charactervector, isovercomplete, isCSS, relativedistance, logicals,
    logicaloperators, barelogicals, bare, logicalsmatrix, gauges, gaugeoperators, gaugesmatrix,
    gaugeoperatorsmatrix, dressed, dressedoperators, dressedlogicals, gaugegroup, gaugegroupmatrix,
    gaugegeneratorsmatrix, gaugegroupgeneratorsmatrix, setsigns!, setlogicals!, setminimumdistance!,
    splitstabilizers, islogical, syndrome, Xsyndrome, Zsyndrome, promotelogicalstogauge!, swapXZlogicals!,
    swapXZgaugeoperators!, allstabilizers, elements, printallstabilizers, printallelements,
    augment, expurgate, fixgauge, setXstabilizers, setZstabilizers, setstabilizers,
    setZstabilizers!, setdistancelowerbound!, permutecode!, permutecode, setstabilizers!,
    setXstabilizers!, standardformA, standardformA1, standardformA2, standardformB, standardformC1,
    standardformC2, standardformD, standardformE, logicalsstandardform, promotegaugestological!,
    promotegaugestological, promotelogicalstogauge

#############################
      # stabilizercode.jl
#############################

include("stabilizercode.jl")
export StabilizerCodeCSS, CSSCode, StabilizerCode

#############################
      # graphstate.jl
#############################

include("graphstate.jl")
export ClusterState, GraphState

#############################
  # miscknownquantumcodes.jl
#############################

include("miscknownquantumcodes.jl")
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
export polynomial, type, CWEtoHWE, weightenumerator, MacWilliamsIdentity,
    weightdistribution, weightplot, support, minimumdistance, weightplotCSSX,
    weightplotCSSZ, weightplotCSS, minimumdistanceXZ, minimumdistanceX,
    minimumdistanceZ, ispure, Sternsattack, Graycodemindist, minimumwords,
    wordsofweight

#############################
   # quantumproductcodes.jl
#############################

include("quantumproductcodes.jl")
export HypergraphProductCode, GeneralizedShorCode, BaconCasaccinoConstruction,
    HyperBicycleCodeCSS, HyperBicycleCode, GeneralizedBicycleCode,
    GeneralizedHypergraphProductCodeMatrices, LiftedGeneralizedHypergraphProductCode,
    QuasiCyclicLiftedProductCodeMatrices, QuasiCyclicLiftedProductCode,
    BiasTailoredQuasiCyclicLiftedProductCodeMatrices, BiasTailoredQuasiCyclicLiftedProductCode

#############################
        # tilings.jl
#############################

include("tilings.jl")
export ReflectionGroup, trianglegroup, rsgroup, tetrahedrongroup, qrsgroup,
    startetrahedrongroup, cycletetrahedrongroup, normalsubgroups, fixedpointfree,
    orientable, kcolorable, cosetintersection

#############################
        # Tanner.jl
#############################

include("Tanner.jl")
export Tannergraphplot, Tannergraph, Tannercode

end
