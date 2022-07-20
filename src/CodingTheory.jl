module CodingTheory

# change environment variable so that banner doesn't print
ENV["NEMO_PRINT_BANNER"] = "false"

using AbstractAlgebra
using Nemo

import LinearAlgebra: tr
import AbstractAlgebra: quo, VectorSpace
import Nemo: isprime, factor, transpose, order, polynomial
import Base: show, length, in, zeros, ⊆, /, *, ==, ∩, +, -

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
      # miscknowncodes.jl
#############################

include("miscknowncodes.jl")

#############################
      # quantumcode.jl
#############################

abstract type AbstractAdditiveCode <: AbstractCode end
abstract type AbstractStabilizerCode <: AbstractAdditiveCode end
abstract type AbstractCSSCode <: AbstractStabilizerCode end

include("quantumcode.jl")

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
        # Exports
#############################

#############################
         # utils.jl
#############################

export kroneckerproduct, Hammingweight, weight, wt, Hammingdistance, distance,
    dist, tr, expandmatrix, symplecticinnerproduct, aresymplecticorthogonal,
    Hermitianinnerproduct, Hermitianconjugatematrix, FpmattoJulia, istriorthogonal,
    printstringarray, printchararray, printsymplecticarray, pseudoinverse,
    quadratictosymplectic, symplectictoquadratic
    #, _processstrings,
    #_Paulistringtosymplectic, _removeempty

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
    words, codewords, elements, isMDS
    # _standardform,

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
    CyclicCode, BCHCode, ReedSolomonCode, complement, ==, ∩, +

#############################
     # miscknowncodes.jl
#############################

export RepetitionCode, Hexacode, HammingCode, TetraCode, SimplexCode,
    GolayCode, ExtendedGolayCode

#############################
      # quantumcode.jl
#############################

export AbstractAdditiveCode, AbstractStabilizerCode, AbstractCSSCode

export field, quadraticfield, length, numqubits, dimension, cardinality,
    rate, signs, Xsigns, Zsigns, stabilizers, symplecticstabilizers,
    Xstabilizers, Zstabilizers, numXstabs, numZstabs, normalizermatrix,
    charactervector, relativedistance, splitstabilizers, isCSS, CSSCode,
    QuantumCode, logicalspace, setlogicals!, changesigns!, Xsyndrome, Zsyndrome,
    syndrome, allstabilizers, elements

#############################
  # miscknownquantumcodes.jl
#############################

export fivequbitcode, Q513, Steanecode, Q713, _Steanecodetrellis, Shorcode, Q913,
    Q412, Q422, Q511, Q823, Q15RM, Q1513, Q1573, triangularsurfacecode,
    rotatedsurfacecode, XZZXsurfacecode, tricolorcode488, tricolorcode666

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

export weightenumeratorC, weightenumerator, weightdistribution, minimumdistance,
    Pauliweightenumerator, Pauliweightenumerator, PWEtoHWE, PWEtoXWE, PWEtoZWE,
    HammingweightenumeratorQ, Hammingweightenumerator, weightenumerator,
    weightdistribution, CWEtoHWE, support, polynomial, type, MacWilliamsIdentity
    # _weightenumeratorBF

end
