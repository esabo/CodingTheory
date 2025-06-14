# Copyright (c) 2023 - 2024 Eric Sabo
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

#############################
# general functions
#############################

"""
$TYPEDSIGNATURES

Return a plot of the EXIT chart for the ensemble given the channel up to a numerical tolerance of `tol`.

# Note
- Run `using Makie` to activate this extension.
"""
function CodingTheory.EXIT_chart_plot(
    E::LDPCEnsemble,
    Ch::AbstractClassicalNoiseChannel;
    tol::Float64 = 1e-9,
)

    @assert isa(Ch, BinaryErasureChannel) "Only BEC is implemented so far"
    x = 0:0.01:1
    c = 1 .- [CodingTheory._poly_eval(x, E.ρ) for x in 1 .- x]
    v = Ch.param .* [CodingTheory._poly_eval(x, E.λ) for x in x]
    Ch ∈ keys(E.density_evo) || CodingTheory._density_evolution!(E, Ch)
    evo_x, evo_y = E.density_evo[Ch]
    ind = findfirst(evo_x .<= tol)
    title = if isnothing(ind)
        "Failed to converge after $(length(evo_x) - 1) iterations (tol = $tol)"
    else
        "Converged after $(ind - 1) iterations (tol = $tol), \$\\varepsilon = $(Ch.param)\$"
    end

    fig = Figure()
    ax = Axis(fig[1, 1], limits = ((0, min(1, Ch.param * 1.2)), (0, 1)), title = title)
    lines!(x, c, label = L"c(x)")
    lines!(v, x, label = L"v_{\varepsilon}^{-1}(x)")
    stairs!(evo_x[1:ind], evo_y[1:ind], linestyle = :dash, linecolor = :black)
    axislegend(position = :rb)
    display(fig)
    return fig

    # p = plot(x, c, label = "\$c(x)\$")
    # plot!(p, v, x, label = "\$v_{\\varepsilon}^{-1}(x)\$")
    # plot!(p, evo_x[1:ind], evo_y[1:ind], seriestype = :steppre, linestyle = :dash,
    #     linecolor = :black, label = false)
    # plot!(p, legend = :bottomright, xlims = (0, min(1, Ch.param * 1.2)), ylims = (0, 1),
    #     title = title)
    # return p
end
