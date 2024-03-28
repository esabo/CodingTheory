module MakieExt

using CodingTheory, Makie, NetworkLayout, CairoMakie, GraphMakie, GLMakie, WGLMakie, GraphPlot, Graphs
import Graphs as Grphs
# import CairoMakie: save

include("Classical/Tanner.jl")
include("Classical/weight_dist.jl")
include("LDPC/codes.jl")
include("LDPC/analysis.jl")
include("Quantum/weight_dist.jl")

end
