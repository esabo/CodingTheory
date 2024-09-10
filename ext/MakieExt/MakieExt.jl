module MakieExt

import CodingTheory

import CodingTheory: Tanner_graph_plot, weight_plot, EXIT_chart_plot, degree_distributions_plot, count_short_cycles_plot, count_elementary_cycles_plot, ACE_spectrum_plot, computation_graph, weight_plot_CSS_X, weight_plot_CSS_Z, weight_plot_CSS

import Makie

import NetworkLayout

import CairoMakie

import GraphMakie

import GLMakie

import WGLMakie

import GraphPlot

import Oscar

import Graphs as Grphs
# import CairoMakie: save

include("Classical/Tanner.jl")
include("Classical/weight_dist.jl")
include("LDPC/codes.jl")
include("LDPC/analysis.jl")
include("Quantum/weight_dist.jl")

end
