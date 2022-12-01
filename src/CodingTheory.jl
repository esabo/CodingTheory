# Copyright (c) 2021, 2022 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

module CodingTheory

# change environment variable so that banner doesn't print
ENV["NEMO_PRINT_BANNER"] = "false"

# using AbstractAlgebra
# using Nemo
using Oscar
using CairoMakie, Graphs
# using Plots
using JLD2
using Combinatorics
using .Threads
using GAP #, GAP_jll
using LinearAlgebra # need to remove? want only mul!
using SparseArrays

import LinearAlgebra: tr, Adjoint
import AbstractAlgebra: quo, VectorSpace
import Nemo: isprime, factor, transpose, order, polynomial, nrows, ncols, degree, isisomorphic, lift
import Base: circshift, reverse, iseven, show, length, in, zeros, ⊆, /, *, ==, ∩, +, -
import CairoMakie: save
import Combinatorics: powerset

# don't want this here
# GAP.Packages.load("LINS");

#############################
        # utils.jl
#############################

include("utils.jl")

#############################
        # cyclotomic.jl
#############################

include("cyclotomic.jl")

#############################
       # linearcode.jl
#############################

abstract type AbstractCode end
abstract type AbstractLinearCode <: AbstractCode end

include("linearcode.jl")

#############################
         # LDPC.jl
#############################

abstract type AbstractLDPCCode <: AbstractLinearCode end

include("LDPC.jl")

#############################
        # LDPCalgs.jl
#############################

include("LDPCalgs.jl")

#############################
    # MatrixProductCode.jl
#############################

abstract type AbstractMatrixProductCode <: AbstractLinearCode end

include("MatrixProductCode.jl")

#############################
      # ReedMuller.jl
#############################

abstract type AbstractReedMullerCode <: AbstractLinearCode end

include("ReedMuller.jl")

#############################
      # cycliccode.jl
#############################

abstract type AbstractCyclicCode <: AbstractLinearCode end
abstract type AbstractBCHCode <: AbstractCyclicCode end
abstract type AbstractReedSolomonCode <: AbstractBCHCode end

include("cycliccode.jl")

#############################
    # quasicycliccode.jl
#############################

abstract type AbstractQuasiCyclicCode <: AbstractLinearCode end

include("quasicycliccode.jl")

#############################
      # miscknowncodes.jl
#############################

include("miscknowncodes.jl")

#############################
 # GeneralizedReedSolomon.jl
#############################

abstract type AbstractGeneralizedReedSolomonCode <: AbstractLinearCode end

include("GeneralizedReedSolomon.jl")

#############################
      # quantumcode.jl
#############################

abstract type AbstractAdditiveCode <: AbstractCode end
abstract type AbstractStabilizerCode <: AbstractAdditiveCode end
abstract type AbstractCSSCode <: AbstractStabilizerCode end

include("quantumcode.jl")

#############################
      # graphstate.jl
#############################

abstract type AbstractGraphState <: AbstractStabilizerCode end
abstract type AbstractGraphStateCSS <: AbstractCSSCode end

include("graphstate.jl")

#############################
  # miscknownquantumcodes.jl
#############################

include("miscknownquantumcodes.jl")

#############################
        # trellis.jl
#############################

include("trellis.jl")

#############################
      # weight_dist.jl
#############################

include("weight_dist.jl")

#############################
   # quantumproductcodes.jl
#############################

abstract type AbstractHypergraphProductCode <: AbstractCSSCode end

include("quantumproductcodes.jl")

#############################
        # tilings.jl
#############################

include("tilings.jl")

#############################
        # Tanner.jl
#############################

include("Tanner.jl")

#############################
        # Exports
#############################

#############################
         # utils.jl
#############################

export kroneckerproduct, Hammingweight, weight, wt, Hammingdistance, distance,
    dist, tr, expandmatrix, symplecticinnerproduct, aresymplecticorthogonal,
    Hermitianinnerproduct, Hermitianconjugatematrix, FpmattoJulia, istriorthogonal,
    printstringarray, printchararray, printsymplecticarray, pseudoinverse,
    quadratictosymplectic, symplectictoquadratic, _removeempty, quadraticresidues,
    digitstoint, isbasis, primitivebasis, #polynomialbasis, monomialbasis,
    normalbasis, dualbasis, complementarybasis, verifydualbasis,
    verifycomplementarybasis, isequivalentbasis, isselfdualbasis,
    isprimitivebasis, isnormalbasis, isextension, polytocircmatrix#,
    #circshift, lift
    #, _processstrings,
    #_Paulistringtosymplectic,

#############################
       # cyclotomic.jl
#############################

export ord, cyclotomiccoset, allcyclotomiccosets, complementqcosets,
    qcosetpairings, qcosetpairings, qcosettable, dualqcosets

#############################
       # linearcode.jl
#############################

export AbstractCode, AbstractLinearCode

export WeightEnumerator, LinearCode, field, length, dimension, cardinality, rate,
    relativedistance, generatormatrix, originalgeneratormatrix, paritycheckmatrix,
    originalparitycheckmatrix, genus, Singletonbound, numbercorrectableerrors,
    encode, syndrome, in, ⊆, ⊂, issubcode, codecomplement, quo, quotient, /, dual,
    Hermitiandual, isequivalent, isselfdual, isselforthogonal, isweaklyselfdual, ⊕,
    directsum, ⊗, kron, tensorproduct, directproduct, productcode, extend, puncture,
    expurgate, augment, shorten, lengthen, uuplusv, Plotkinconstruction, subcode,
    juxtaposition, constructionX, constructionX3, upluswvpluswuplusvplusw,
    entrywiseproductcode, *, Schurproductcode, Hadamardproductcode,
    componentwiseproductcode, VectorSpace, setminimumdistance!,
    expandedcode, subfieldsubcode, tracecode, evensubcode, permutecode,
    words, codewords, elements, isMDS, iseven, isdoublyeven, istriplyeven
    # _standardform,

export dual1, dual2, dual3

#############################
         # LDPC.jl
#############################

export AbstractLDPCCode

export variabledegreedistribution, checkdegreedistribution,
    degreedistributions, columnbound, rowbound, bounds, density, isregular,
    LDPCCode, degreedistributionssplot, variabledegreepolynomial,
    checkdegreepolynomial, columnrowbounds

#############################
        # LDPCalgs.jl
#############################

# export Sternsattack

#############################
    # MatrixProductCode.jl
#############################

export AbstractMatrixProductCode

export MatrixProductCode

#############################
       # ReedMuller.jl
#############################

export AbstractReedMullerCode

export order, RMr, RMm, ReedMullergeneratormatrix, ReedMullerCode

#############################
       # cycliccode.jl
#############################

export AbstractCyclicCode, AbstractBCHCode, AbstractReedSolomonCode

export definingset, splittingfield, polynomialring, primitiveroot, offset,
    designdistance, qcosets, qcosetsreps, generatorpolynomial, paritycheckpolynomial,
    idempotent, isprimitive, isnarrowsense, isreversible, finddelta, dualdefiningset,
    CyclicCode, BCHCode, ReedSolomonCode, complement, ==, ∩, +, QuadraticResidueCode

#############################
    # quasicycliccode.jl
#############################

export AbstractQuasiCyclicCode

export weightmatrix, basematrix, protographmatrix, QuasiCyclicCode, index,
    expansionfactor, type, polynomialmatrix, polynomialmatrixtype,
    noncirculantgeneratormatrix, noncirculantparitycheckmatrix, generators,
    circulants, issinglegenerator

#############################
     # miscknowncodes.jl
#############################

export RepetitionCode, Hexacode, HammingCode, TetraCode, SimplexCode,
    GolayCode, ExtendedGolayCode

#############################
 # GeneralizedReedSolomon.jl
#############################

export AbstractGeneralizedReedSolomonCode

export GeneralizedReedSolomonCode, scalars, dualscalars, evaluationpoints

#############################
      # quantumcode.jl
#############################

export AbstractAdditiveCode, AbstractStabilizerCode, AbstractCSSCode

export field, quadraticfield, length, numqubits, dimension, cardinality,
    rate, signs, Xsigns, Zsigns, stabilizers, symplecticstabilizers,
    Xstabilizers, Zstabilizers, numXstabs, numZstabs, normalizermatrix,
    charactervector, relativedistance, splitstabilizers, isCSS, CSSCode,
    QuantumCode, logicalspace, setlogicals!, changesigns!, Xsyndrome, Zsyndrome,
    syndrome, allstabilizers, elements, logicals, logicalsmatrix, isisomorphic

#############################
      # graphstate.jl
#############################

export AbstractGraphState, AbstractGraphStateCSS

export ClusterState, graphstate

#############################
  # miscknownquantumcodes.jl
#############################

export FiveQubitCode, Q513, SteaneCode, Q713, _SteaneCodeTrellis, ShorCode, Q913,
    Q412, Q422, Q511, Q823, Q15RM, Q1513, Q1573, TriangularSurfaceCode,
    RotatedSurfaceCode, XZZXSurfaceCode, TriangularColorCode488, TriangularColorCode666,
    ToricCode, PlanarSurfaceCode, XYSurfaceCode, XYZ2Code

#############################
        # trellis.jl
#############################

export Trellis, vertices, edges, isisomorphic, isequal, loadbalancedecode,
    trellisorientedformlinear,trellisprofiles, syndrometrellis,
    trellisorientedformadditive, optimalsectionalizationQ, weightQ!,
    shiftandweightQ!, shiftanddecodeQ!, shift!, isshifted

#############################
      # weight_dist.jl
#############################

export polynomial, type, CWEtoHWE, weightenumerator, MacWilliamsIdentity,
    weightdistribution, weightplot, support, minimumdistance, weightplotCSSX,
    weightplotCSSZ, weightplotCSS, minimumdistanceXZ, minimumdistanceX,
    minimumdistanceZ, ispure, Sternsattack, Graycodemindist

#############################
   # quantumproductcodes.jl
#############################

export HypergraphProductCode, GeneralizedShorCode, BaconCasaccinoConstruction,
    HyperBicycleCodeCSS, HyperBicycleCode, GeneralizedBicycleCode,
    BicycleCode, GeneralizedHypergraphProductCode, LiftedGeneralizedHypergraphProductCode,
    QuasiCyclicLiftedProductCode, LiftedQuasiCyclicLiftedProductCode,
    BiasTailoredQuasiCyclicLiftedProductCode, LiftedBiasTailoredQuasiCyclicLiftedProductCode

#############################
        # tilings.jl
#############################

export ReflectionGroup, trianglegroup, rsgroup, tetrahedrongroup, qrsgroup,
    startetrahedrongroup, cycletetrahedrongroup, normalsubgroups, fixedpointfree,
    orientable, kcolorable, cosetintersection

#############################
        # Tanner.jl
#############################

export Tannergraph, Tannercode

end