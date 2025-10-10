# Copyright (c) 2025 Eric Sabo, David Marquis
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
      # abstract types
#############################

abstract type AbstractLattice end

#############################
      # concrete types
#############################

struct Lattice <: AbstractLattice
    B::MatrixElem{<:RingElem}
end
