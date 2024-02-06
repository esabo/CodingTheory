# CodingTheory

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://esabo.github.io/CodingTheory/dev/)
[![Build Status](https://github.com/esabo/CodingTheory/actions/workflows/Tests.yml/badge.svg?branch=master)](https://github.com/esabo/CodingTheory/actions/workflows/Tests.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/esabo/CodingTheory/branch/master/graph/badge.svg)](https://codecov.io/gh/esabo/CodingTheory)

A coding theory library for Julia.

The goal of this package is to develop a classical and quantum error-correcting codes package in as much native Julia as possible. The library is built around the Oscar.jl framework, and many thanks to Tommy Hofmann of these packages for helping this repo get off the ground. Anyone is welcome to contribute, although the final form of any accepted code may be standardized to maintain intra-package consistency.

At the moment, all functions work as intended for test cases but have not been unit tested thoroughly enough to guarantee accuracy and error free usage. All results from this library should be mentally checked and any bugs reported (or fixed and pushed). This is particularly true for the quantum part where algorithms become increasingly more complicated.

The generation of hyperbolic tilings (in tilings.jl) requires the [LINS package](https://github.com/FriedrichRober/LINS). Please use the [temporary fork](https://github.com/esabo/LINS.git) at the address below to fix a compatibility issue with the Sonata package. The probabilistic quantum minimum distance algorithm [QDistRnd](https://github.com/QEC-pages/QDistRnd.git) requires the corresponding GAP package. To install these, run 
    GAP.Packages.install(URL)
    GAP.Packages.install("https://github.com/QEC-pages/QDistRnd.git")
    GAP.Packages.install("https://github.com/esabo/LINS.git")
and then
    GAP.Packages.load("QDistRnd");
    GAP.Packages.load("LINS");
These are not automatically installed and loaded.

Quantum minimum distance functions are currently disabled due to a change in the underlying structs. DistRandCSS is still available via the built-in GAP interface.

Parts of the library are multi-threaded and benefit greatly from the use of multiple cores.

A growing list of examples and tutorials are provided in the documentation (which is only slightly out of date).