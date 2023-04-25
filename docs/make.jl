using Documenter, CodingTheory

# change to trigger bot

# root = "../",
# 	source = "src",
# 	build = "build",
Documenter.makedocs(
	clean = true,
	doctest = false,
	modules = Module[CodingTheory],
	repo = "",
	highlightsig = true,
	sitename = "Coding Theory Documentation",
	expandfirst = [],
    pages = ["Table Of Contents" => "index.md",
		     "Linear Codes" => "linearcode.md",
		     "Cyclic Codes" => "cycliccode.md",
			 "Generalized Reed-Solomon Codes" => "GeneralizedReedSolomon.md",
			 "Quasi-Cyclic Codes" => "quasicyclic.md",
		     "Reed-Muller Codes" => "ReedMuller.md",
		     "Miscellaneous Known Codes" => "miscknowncodes.md",
			 "LDPC Codes" => "LDPC.md",
			 "Tilings Of Surfaces" => "tilings.md",
		     "Stabilizer Codes" => "quantumcode.md",
		     "Miscellaneous Known Stabilizer Codes" => "miscknownquantumcodes.md",
			 "Product Codes" => "productcodes.md",
		     "Trellises" => "trellis.md",
		     "Weight Enumerators, Distributions, And Minimum Distances" => "weight_dist.md",
		     "Utilies" => "utils.md",
             "Index" => "theindex.md"]
)

deploydocs(repo = "github.com/esabo/CodingTheory.git")
