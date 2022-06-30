module CodingTheory


####################
# LOAD PACKAGES
####################

using AbstractAlgebra

# change environment variable so that banner doesn't print
ENV["NEMO_PRINT_BANNER"] = "false"
using Nemo

import Base: show, length, in, ⊆, /, *, ==, ∩, +
import AbstractAlgebra: quo, VectorSpace

# I want to remove this dependence
using Primes
import Primes: factor

import LinearAlgebra: tr

using Reexport

# using SymPy
# using Plots


#####################
# TYPE STRUCTURE
#####################

abstract type AbstractCode end
abstract type AbstractLinearCode <: AbstractCode end
abstract type AbstractCyclicCode <: AbstractLinearCode end
abstract type AbstractBCHCode <: AbstractCyclicCode end
abstract type AbstractReedSolomonCode <: AbstractBCHCode end
abstract type AbstractReedMullerCode <: AbstractLinearCode end

# can one day do classical additive codes
abstract type AbstractAdditiveCode <: AbstractCode end
abstract type AbstractStabilizerCode <: AbstractAdditiveCode end
abstract type AbstractCSSCode <: AbstractStabilizerCode end
# can also build a AbstractQuantumLinearCode if the additive code is also linear
# will probably have to tweak all of these later

export AbstractCode, AbstractLinearCode, AbstractCyclicCode, AbstractBCHCode, AbstractReedSolomonCode, AbstractReedMullerCode, AbstractAdditiveCode, AbstractStabilizerCode, AbstractCSSCode


######################
# LOAD FILES
######################
# module Utils
include("utils.jl")
# export asdf
# end
# @reexport using CodingTheory.Utils

include("cyclotomic.jl")
include("linearcode.jl")
include("ReedMuller.jl")
include("cycliccode.jl")
include("miscknowncodes.jl")
include("quantumcode.jl")
include("miscknownquantumcodes.jl")
include("trellis.jl")
include("weight_dist.jl")

# TODO: should be in a better place
export fq_nmod_mat

#######################
# EXPORTS
#######################

#############################
       # linearcode.jl
#############################

export LinearCode, field, length, dimension, cardinality, rate, setminimumdistance,
relativedistance, generatormatrix, originalgeneratormatrix, paritycheckmatrix,
originalparitycheckmatrix, genus, Singletonbound, numbercorrectableerrors,
encode, syndrome, in, ⊆, ⊂, issubcode, codecomplement, quo, quotient, /, dual,
Hermitiandual, isequivalent, isselfdual, isselforthogonal, isweaklyselfdual, ⊕,
directsum, ⊗, kron, tensorproduct, directproduct, productcode, extend, puncture,
expurgate, augment, shorten, lengthen, uuplusv, Plotkinconstruction, subcode,
juxtaposition, constructionX, constructionX3, upluswvpluswuplusvplusw,
expandedcode, entrywiseproductcode, *, Schurproductcode, Hadamardproductcode,
componentwiseproductcode, VectorSpace

#############################
       # cycliccode.jl
#############################

export definingset, splittingfield, polynomialring, primitiveroot, offset,
designdistance, qcosets, qcosetsreps, generatorpolynomial, paritycheckpolynomial,
idempotent, isprimitive, isnarrowsense, isreversible, finddelta, dualdefiningset,
CyclicCode, BCHCode, ReedSolomonCode, complement, ==, ∩, +

#############################
       # cyclotomic.jl
#############################

export ord, cyclotomiccoset, allcyclotomiccosets, complementqcosets,
qcosetpairings, qcosetpairings, qcosettable, dualqcosets

#############################
     # miscknowncodes.jl
#############################

export repetitioncode

#############################
 # miscknownquantumcodes.jl
#############################

export fivequbitcode, Q513, Steanecode, Q713, _Steanecodetrellis, Shorcode, Q913,
Q412, Q422, Q511, Q823, Q15RM, Q1513, Q1573, triangularsurfacecode,
rotatedsurfacecode, XZZXsurfacecode, tricolorcode488, tricolorcode666

#############################
      # ReedMuller.jl
#############################

export order, RMr, RMm, ReedMullergeneratormatrix, ReedMullerCode

#############################
        # trellis.jl
#############################

export vertices, edges, isisomorphic, isequal, loadbalancedecode,
trellisorientedformC,trellisprofiles, syndrometrellis, trellisorientedformQ,
optimalsectionalizationQ, weightQ!, shiftandweightQ!, shiftanddecodeQ!,
shift!

#############################
     # weight_dist.jl
#############################

export weightenumeratorC, weightenumerator, weightdistribution, minimumdistance,
Pauliweightenumerator, Pauliweightenumerator, PWEtoHWE, PWEtoXWE, PWEtoZWE,
HammingweightenumeratorQ, Hammingweightenumerator, weightenumerator,
weightdistribution

#############################
         # utils.jl
#############################

export kroneckerproduct, Hammingweight, weight, wt, Hammingdistance, distance,
dist, tr, expandmatrix, symplecticinnerproduct, aresymplecticorthogonal,
Hermitianinnerproduct, Hermitianconjugatematrix, FpmattoJulia, istriorthogonal,
printstringarray, printchararray, printsymplecticarray, pseudoinverse,
quadratictosymplectic, symplectictoquadratic


end # module
