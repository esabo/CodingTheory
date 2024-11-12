# CodingTheory

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://esabo.github.io/CodingTheory/dev/)
[![Build Status](https://github.com/esabo/CodingTheory/actions/workflows/Tests.yml/badge.svg?branch=master)](https://github.com/esabo/CodingTheory/actions/workflows/Tests.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/esabo/CodingTheory/branch/master/graph/badge.svg)](https://codecov.io/gh/esabo/CodingTheory)

A coding theory library for Julia.

The goal of this package is to develop a classical and quantum error-correcting codes package in as much native Julia as possible. The library is built around the Oscar.jl framework, and many thanks to Tommy Hofmann of these packages for helping this repo get off the ground. Anyone is welcome to contribute, although the final form of any accepted code may be standardized to maintain intra-package consistency.

At the moment, all functions work as intended for test cases but have not been unit tested thoroughly enough to guarantee 100% accuracy and error free usage. All results from this library should be mentally checked and any bugs reported (or fixed and pushed).

Parts of the library are multi-threaded and benefit greatly from the use of multiple cores.

The minimum distance functions are currently being rewritten and will reappear soon. Depending on what one is looking for, current functions may be sufficient. Feel free to reach out on the [Slack channel](https://join.slack.com/t/juliacodingtheory/shared_invite/zt-2u8n5h5wm-QqnXl2NZqRvTmGGEPumbqQ).

Improved documentation is currently a major to-do. Again, feel free to ask questions on Slack.
