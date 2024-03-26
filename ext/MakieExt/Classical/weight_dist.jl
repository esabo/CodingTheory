# Copyright (c) 2022 - 2024 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
    # Weight Enumerators
#############################

"""
    weight_plot(C::AbstractLinearCode; alg::Symbol = :auto)

Return a bar graph of the weight distribution of `C`.

# Note
- Run `using Makie` to activate this extension.
"""
function weight_plot(C::AbstractLinearCode; alg::Symbol = :auto)
    wt_dist = weight_distribution(C, alg = alg, compact = true)
    x_ticks = findall(x -> x > 0, vec(wt_dist)) .- 1
    y_ticks = [wt_dist[i] for i in 1:length(wt_dist) if !iszero(wt_dist[i])]
    ismissing(C.d) ? (title = "Weight Distribution - [$(C.n), $(C.k)]";) :
        title = "Weight Distribution - [$(C.n), $(C.k), $(C.d)]"

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Weight", ylabel = "Number of Terms", title = title)
    barplot!(ax, 0:C.n, wt_dist', bar_width = 1, xticks = x_ticks, yticks = y_ticks)
    # fig = bar(0:C.n, wt_dist', bar_width = 1, xticks = x_ticks, yticks = y_ticks,
    #     legend = false, xlabel = "Weight", ylabel = "Number of Terms", title = title)
    display(fig)
    return fig
end
