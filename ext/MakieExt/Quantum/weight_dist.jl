# Copyright (c) 2022 - 2024 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
    # Weight Enumerators
#############################

# """
# $TYPEDSIGNATURES

# Return a bar graph of the weight distribution related to `S`.

# If `type` is `:stabilizer`, the weight distribution of the stabilizers are computed.
# If `type` is `:normalizer`, the weight distrbution of the normalizer of the stabilizers
# are computed. If `type` is `:quotient`, the weight distrbution of the normalizer mod the
# stabilizers is computed.

# # Note
# - Run `using Makie` to activate this extension.
# """
function CodingTheory.weight_plot(S::AbstractStabilizerCode; alg::Symbol = :auto,
    type::Symbol = :stabilizer)

    type ∈ (:stabilizer, :normalizer, :quotient) || throw(ArgumentError("Unknown value $type for parameter type."))

    wt_dist = weight_distribution(S, alg = alg, type = type, compact = false)
    x_ticks = findall(x -> x > 0, vec(wt_dist)) .- 1
    y_ticks = [wt_dist[i] for i in 1:length(wt_dist) if !iszero(wt_dist[i])]
    if type == :stabilizer
        title_str = "Stabilizer Weight Distribution"
    elseif type == :normalizer
        title_str = "Normalizer Weight Distribution"
    else
        title_str = "Quotient Weight Distribution"
    end
    ismissing(S.d) ? (title = "$title_str - [$(S.n), $(S.k)]";) :
        title = "$title_str - [$(S.n), $(S.k), $(S.d)]"

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Weight", ylabel = "Number of Terms", title = title)
    barplot!(ax, 0:S.n, wt_dist', bar_width = 1, xticks = x_ticks, yticks = y_ticks)
    # fig = bar(0:S.n, wt_dist', bar_width = 1, xticks = x_ticks, yticks = y_ticks,
    #     legend = false, xlabel = "Weight", ylabel = "Number of Terms", title = title)
    # display(fig)
    return fig
end

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
