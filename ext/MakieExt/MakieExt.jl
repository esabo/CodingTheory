module MakieExt

import CodingTheory, Oscar
import CodingTheory: Tanner_graph_plot, weight_plot, EXIT_chart_plot, degree_distributions_plot, simple_cycle_length_distribution_plot, ACE_spectrum_plot, computation_graph, weight_plot_CSS_X, weight_plot_CSS_Z, weight_plot_CSS, AbstractLinearCode, AbstractLDPCCode,
    LDPCEnsemble, AbstractClassicalNoiseChannel, AbstractStabilizerCode, AbstractStabilizerCodeCSS,
    simple_cycle_length_distribution, simple_cycle_distribution_by_variable_node, ACE_distribution,
    short_cycle_length_distribution, short_cycle_distribution_by_variable_node, girth
# import Makie
using Makie, NetworkLayout, CairoMakie, GraphMakie, GLMakie, WGLMakie#, GraphPlot
import CairoMakie: @L_str #, Figure, Axis, hidespines!, hidedecorations!, lines!, text!
import Graphs as Grphs
using DocStringExtensions
# import CairoMakie: save

include("Classical/Tanner.jl")
include("Classical/weight_dist.jl")
include("LDPC/codes.jl")
include("LDPC/cycles.jl")
include("LDPC/analysis.jl")
include("Quantum/weight_dist.jl")

end
