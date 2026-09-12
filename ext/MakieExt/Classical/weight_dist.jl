# Copyright (c) 2022 - 2026 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
    # Weight Enumerators
#############################

"""
$TYPEDSIGNATURES

Return a bar graph of the weight distribution of `C`.

# Note
- Run `using Makie` to activate this extension.
"""
function CodingTheory.weight_plot(C::AbstractLinearCode; verbose::Bool=false)
    # Fetch the dense n + 1 array
    wt_dist = CodingTheory.weight_distribution_array(C; verbose=verbose)
    
    # Extract non-zero weights and their counts
    x_ticks = Int[]
    y_ticks = BigInt[]
    
    for (idx, count) in enumerate(wt_dist)
        if !iszero(count)
            push!(x_ticks, idx - 1) # Shift back to 0-based Hamming weight
            push!(y_ticks, count)
        end
    end

    ismissing(C.d) ?
        (title = "Weight Distribution - [$(C.n), $(C.k)]";) :
        title = "Weight Distribution - [$(C.n), $(C.k), $(C.d)]"

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Weight", ylabel = "Number of Terms", title = title)
    
    # Convert ticks to standard numeric types Makie expects
    barplot!(ax, x_ticks, Float64.(y_ticks), bar_width = 1)
    
    display(fig)
    return fig
end