# CodingTheory

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://esabo.github.io/CodingTheory/dev/)
[![Build Status](https://github.com/esabo/CodingTheory/actions/workflows/Tests.yml/badge.svg?branch=dev)](https://github.com/esabo/CodingTheory/actions/workflows/Tests.yml?query=branch%3Adev)
[![Coverage](https://codecov.io/gh/esabo/CodingTheory/graph/badge.svg?branch=dev)](https://codecov.io/gh/esabo/CodingTheory/tree/dev)

A classical, LDPC, and quantum coding theory library for Julia.

The library uses [Oscar.jl](https://www.oscar-system.org/) for exact
finite-field and polynomial arithmetic and native Julia data structures for
performance-sensitive sparse and iterative algorithms.

Install the development version from Julia's package prompt:

```julia
] add https://github.com/esabo/CodingTheory
```

Then:

```julia
using Oscar
using CodingTheory
```

See the [development documentation](https://esabo.github.io/CodingTheory/dev/)
for tutorials and API references. Exact minimum-distance computations can be
exponential; the documentation explains solver selection and the distinction
between certified bounds and witnessed or heuristic results.

Parts of the library are multi-threaded and benefit from multiple Julia
threads. Questions and development discussion are welcome on the
[CodingTheory Slack channel](https://join.slack.com/t/juliacodingtheory/shared_invite/zt-2u8n5h5wm-QqnXl2NZqRvTmGGEPumbqQ).
