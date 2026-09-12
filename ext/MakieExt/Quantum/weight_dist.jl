# Copyright (c) 2022 - 2024 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
$TYPEDSIGNATURES

Return a bar graph of the weight distribution of the `X` stabilizers.

# Note
- Run `using Makie` to activate this extension.
"""
function CodingTheory.weight_plot_CSS_X(S::AbstractStabilizerCodeCSS; alg::Symbol = :auto)
    C = LinearCode(S.X_stabs)
    wt_dist = weight_distribution(C, alg = alg, compact = false)
    x_ticks = findall(x -> x > 0, vec(wt_dist)) .- 1
    y_ticks = [wt_dist[i] for i in 1:length(wt_dist) if !iszero(wt_dist[i])]
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Weight", ylabel = "Number of Terms",
        title = "X-Weight Distribution")
    barplot!(ax, 0:C.n, wt_dist', bar_width = 1, xticks = x_ticks, yticks = y_ticks)
    # fig = bar(0:C.n, wt_dist', bar_width = 1, xticks = x_ticks, yticks = y_ticks,
    #     legend = false, xlabel = "Weight", ylabel = "Number of Terms",
    #     title = "X-Weight Distribution")
    display(fig)
    return fig
end

"""
$TYPEDSIGNATURES

Return a bar graph of the weight distribution of the `Z` stabilizers.

# Note
- Run `using Makie` to activate this extension.
"""
function CodingTheory.weight_plot_CSS_Z(S::AbstractStabilizerCodeCSS; alg::Symbol = :auto)
    C = LinearCode(S.Z_stabs)
    wt_dist = weight_distribution(C, alg = alg, compact = false)
    x_ticks = findall(x -> x > 0, vec(wt_dist)) .- 1
    y_ticks = [wt_dist[i] for i in 1:length(wt_dist) if !iszero(wt_dist[i])]
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Weight", ylabel = "Number of Terms",
        title = "Z-Weight Distribution")
    barplot!(ax, 0:C.n, wt_dist', bar_width = 1, xticks = x_ticks, yticks = y_ticks)
    # fig = bar(0:C.n, wt_dist', bar_width = 1, xticks = x_ticks, yticks = y_ticks,
    #     legend = false, xlabel = "Weight", ylabel = "Number of Terms",
    #     title = "Z-Weight Distribution")
    display(fig)
    return fig
end

"""
$TYPEDSIGNATURES

Return bar graphs of the weight distribution of both the `X` and 'Z' stabilizers, separately.

# Note
- Run `using Makie` to activate this extension.
"""
function CodingTheory.weight_plot_CSS(S::AbstractStabilizerCodeCSS; alg::Symbol = :auto)
    C = LinearCode(S.X_stabs)
    wt_dist = weight_distribution(C, alg = alg, compact = false)
    x_ticks = findall(x -> x > 0, vec(wt_dist)) .- 1
    y_ticks = [wt_dist[i] for i in 1:length(wt_dist) if !iszero(wt_dist[i])]
    fig = Figure()
    ax1 = Axis(fig[1, 1], xlabel = "Weight", ylabel = "Number of Terms",
        title = "X-Weight Distribution")
    barplot!(ax1, 0:C.n, wt_dist', bar_width = 1, xticks = x_ticks, yticks = y_ticks)
    # f_X = bar(0:C.n, wt_dist', bar_width = 1, xticks = x_ticks, yticks = y_ticks,
    #     legend = false, xlabel = "Weight", ylabel = "Number of Terms",
    #     title = "X-Weight Distribution")

    # okay to overwrite
    C = LinearCode(S.Z_stabs)
    wt_dist = weight_distribution(C, alg = alg, compact = false)
    x_ticks = findall(x -> x > 0, vec(wt_dist)) .- 1
    y_ticks = [wt_dist[i] for i in 1:length(wt_dist) if !iszero(wt_dist[i])]
    ax2 = Axis(fig[1, 2], xlabel = "Weight", ylabel = "Number of Terms",
        title = "Z-Weight Distribution")
    barplot!(ax2, 0:C.n, wt_dist', bar_width = 1, xticks = x_ticks, yticks = y_ticks)
    display(fig)
    return fig
    # f_Z = bar(0:C.n, wt_dist', bar_width = 1, xticks = x_ticks, yticks = y_ticks,
    #     legend = false, xlabel = "Weight", ylabel = "Number of Terms",
    #     title = "Z-Weight Distribution")
    # f = Plots.plot(f_X, f_Z, layout = (1, 2))
    # display(f)
    # return f
end
