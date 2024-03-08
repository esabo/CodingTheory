# Copyright (c) 2023 - 2024 Eric Sabo, Benjamin Ide
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
        # LP Decoders
#############################

# function _init_LP_decoder_LDPC end

# function _LP_decoder_LDPC end

# TODO: docstring and in extension
"""
    LP_decoder_LDPC(H::Union{CTMatrixTypes, AbstractMatrix{<:Number}}, v::Union{CTMatrixTypes, Vector{<:Integer}}, Ch::BinarySymmetricChannel)
    LP_decoder_LDPC(C::AbstractLinearCode, v::Union{CTMatrixTypes, Vector{<:Integer}}, Ch::BinarySymmetricChannel)

Return

# Note
- Run `using JuMP, GLPK` to activate this extension.
"""
function LP_decoder_LDPC end
