# Copyright (c) 2021, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

include("linearcode.jl")

abstract type AbstractHammingCode <: AbstractLinearCode end

function construct_ham_matrix(r::Int, q::Int)
    ncols = Int(floor((q^r - 1) / (q - 1)))
    M = Matrix{Int}(undef, r, ncols)

    for i in 1:ncols
        M[:, i] = reverse(digits(parse(Int, string(i, base = q)), pad = r), dims = 1)
    end

    return M
end
