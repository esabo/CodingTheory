# Copyright (c) 2021, Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

include("linearcode.jl")

abstract type AbstractGolayCode <: AbstractLinearCode end

# both of these are actually the extended Golay codes
function Golay(p::Integer)
    if p == 2
        F, _ = FiniteField(2, 1, "α")
        M = MatrixSpace(F, 12 , 12)
        A = M([0 1 1 1 1 1 1 1 1 1 1 1;
             1 1 1 0 1 1 1 0 0 0 1 0;
             1 1 0 1 1 1 0 0 0 1 0 1;
             1 0 1 1 1 0 0 0 1 0 1 1;
             1 1 1 1 0 0 0 1 0 1 1 0;
             1 1 1 0 0 0 1 0 1 1 0 1;
             1 1 0 0 0 1 0 1 1 0 1 1;
             1 0 0 0 1 0 1 1 0 1 1 1;
             1 0 0 1 0 1 1 0 1 1 1 0;
             1 0 1 0 1 1 0 1 1 1 0 0;
             1 1 0 1 1 0 1 1 1 0 0 0;
             1 0 1 1 0 1 1 1 0 0 0 1])
        G = hcat(M(1), A)
        H = hcat(-transpose(A), M(1))
        return LinearCode(F, 24, 12, 8, G, G, H, H, G, H, missing)
    elseif p == 3
        F, _ = FiniteField(3, 1, "α")
        M = MatrixSpace(F, 6 , 6)
        A = M([0 1 1 1 1 1;
               1 0 1 -1 -1 1;
               1 1 0 1 -1 -1;
               1 -1 1 0 1 -1;
               1 -1 -1 1 0 1;
               1 1 -1 -1 1 0])
        G = hcat(M(1), A)
        H = hcat(-transpose(A), M(1))
        return LinearCode(F, 12, 6, 6, G, G, H, H, G, H, missing)
    else
        error("Golay code not implemented for q = $q.")
    end
end
