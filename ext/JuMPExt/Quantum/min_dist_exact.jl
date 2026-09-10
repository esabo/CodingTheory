# Copyright (c) 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

include("../../shared/css_min_dist_ilp.jl")

function _minimum_distance_css_HiGHS(H, logical_checks; kwargs...)
    return _solve_css_minimum_distance_ilp(
        HiGHS.Optimizer, H, logical_checks;
        verbose_attributes=[
            "output_flag" => true,
            "log_to_console" => true,
        ], kwargs...)
end
