# Copyright (c) 2023 - 2025 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
# simple cycles
#############################

"""
$TYPEDSIGNATURES

Return a bar graph and dictionary of (length, count) pairs for the unique simple cycles up to
length `len` of the Tanner graph of `L`. If `len` is `-1`, then all simple cycles will be
enumerated. An empty figure and dictionary are returned when there are no cycles.

# Note
- Simple cycles do not contain the same vertex twice.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
- Run `using Makie` to activate this extension.
"""
function CodingTheory.simple_cycle_length_distribution_plot(
    L::AbstractLDPCCode;
    len::Int = -1,
)
    dist = simple_cycle_length_distribution(L, len = len)
    x_data = collect(keys(dist))
    y_data = collect(values(dist))

    fig = Figure();
    ax = Axis(
        fig[1, 1],
        xlabel = "Cycle Length",
        ylabel = "Occurrences",
        title = "Simple Cycle Counts",
        xticks = x_data,
        yticks = y_data,
    )
    barplot!(ax, x_data, y_data, strokewidth = 1)
    display(fig)
    return fig
end

"""
$TYPEDSIGNATURES

Return bar graph and a dictionary of (node, count) pairs for the unique simple cycles up to length
`len` of the Tanner graph of `L`. If `len` is `-1`, then all simple cycles will be enumerated. An
empty figure and dictionary are returned when there are no cycles.

# Note
- Simple cycles do not contain the same vertex twice.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
    already cached.
- Run `using Makie` to activate this extension.
"""
function CodingTheory.simple_cycle_distribution_by_variable_node_plot(
    L::AbstractLDPCCode;
    len::Int = -1,
)

    dist = simple_cycle_distribution_by_variable_node(L, len = len)
    collect(keys(dist))
    collect(values(dist))

    fig = Figure();
    ax = Axis(
        fig[1, 1],
        xlabel = "Variable Node Index",
        ylabel = "Occurrences",
        title = "Simple Cycles By Variable Node",
        xticks = x_data,
        yticks = y_data,
    )
    barplot!(ax, x_data, y_data, strokewidth = 1)
    display(fig)
    return fig
end

#############################
# short cycles
#############################

"""
$TYPEDSIGNATURES

Return a bar graph and dictionary of (length, count) pairs for the unique short cycles up to
length `len` of the Tanner graph of `L`. If `len` is `-1`, then all short cycles will be
enumerated. An empty figure and dictionary are returned when there are no cycles.

# Note
- Short cycles are defined to be those with lengths between ``g`` and ``2g - 2``,
  where ``g`` is the girth.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
  already cached.
- Run `using Makie` to activate this extension.
"""
function CodingTheory.short_cycle_length_distribution_plot(
    L::AbstractLDPCCode;
    len::Int = -1,
)
    dist = short_cycle_length_distribution(L, len = len)
    x_data = collect(keys(dist))
    y_data = collect(values(dist))

    fig = Figure()
    ax = Axis(
        fig[1, 1],
        xlabel = "Cycle Length",
        ylabel = "Occurrences",
        title = "Short Cycle Counts",
        xticks = x_data,
        yticks = y_data,
    )
    barplot!(ax, x_data, y_data, strokewidth = 1)
    display(fig)
    return fig
end

"""
$TYPEDSIGNATURES

Return bar graph and a dictionary of (node, count) pairs for the unique short cycles up to length
`len` of the Tanner graph of `L`. If `len` is `-1`, then all short cycles will be enumerated. An
empty figure and dictionary are returned when there are no cycles.

# Note
- Short cycles are defined to be those with lengths between ``g`` and ``2g - 2``,
  where ``g`` is the girth.
- This function calls `enumerate_simple_cycles(L, len = len)`, which could be expensive if not
    already cached.
- Run `using Makie` to activate this extension.
"""
function CodingTheory.short_cycle_distribution_by_variable_node_plot(
    L::AbstractLDPCCode;
    len::Int = -1,
)

    dist = short_cycle_distribution_by_variable_node(L, len = len)
    x_data = collect(keys(dist))
    y_data = collect(values(dist))

    fig = Figure();
    ax = Axis(
        fig[1, 1],
        xlabel = "Variable Node Index",
        ylabel = "Occurrences",
        title = "Short Cycles By Variable Node",
        xticks = x_data,
        yticks = y_data,
    )
    barplot!(ax, x_data, y_data, strokewidth = 1)
    display(fig)
    return fig
end

#############################
# lollipops
#############################


# one for short cycles which cuts off stems at the correct lengths
# one for doing all of them





#############################
# ACE
#############################

"""
$TYPEDSIGNATURES

Return an interactive figure and data for the ACE spectrum of the Tanner graph of `C`.

# Note
- Run `using Makie` to activate this extension.
"""
function CodingTheory.ACE_spectrum_plot(C::AbstractLDPCCode)
    vs_ACEs, _ = ACE_distribution(C, collect(1:C.n))
    grth = girth(C)
    # TODO: remove WGLMakie as a default use and only use for interactive plots
    fig = Figure();
    ax = Axis(fig[1, 1], xlabel = "ACE", ylabel = "Occurrences", title = "ACE Spectrum")
    sg = SliderGrid(
        fig[2, 1],
        (label = "Cycle Length", range = grth:2:(2*grth-2), startvalue = 4),
    )

    x_max = maximum(reduce(vcat, vs_ACEs))
    y_max = 0
    counts = ACE_spectrum(C)
    X_data = Observable(Vector{Int}())
    Y_data = Observable(Vector{Int}())
    barplot!(ax, X_data, Y_data, strokewidth = 1, xticks = X_data, yticks = Y_data)
    for (k, l) in enumerate(grth:2:(2*grth-2))

        ys = collect(values(counts[k]))
        y_max = maximum([y_max; ys])

        on(sg.sliders[1].value) do val
            if to_value(val) == l
                X_data.val = collect(keys(counts[k]))
                Y_data.val = ys
                notify(X_data)
                notify(Y_data)
            end
        end
    end

    GLMakie.limits!(0, x_max + 1, 0, y_max + 2)
    display(fig)
    return fig, counts
end
