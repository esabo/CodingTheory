module MakieExt

using CodingTheory, Makie, NetworkLayout#, GraphPlot
# import CairoMakie: save

include("Classical/Tanner.jl")
include("Classical/weight_dist.jl")
include("LDPC/codes.jl")
include("LDPC/analysis.jl")
include("Quantum/weight_dist.jl")

end
