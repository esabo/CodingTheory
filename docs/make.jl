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
    pages = ["Introduction" => "index.md",
		"Tutorials" => [
			"Tutorials/Linear Codes.md",
			"Tutorials/Cyclic Codes.md",
			"Tutorials/Quantum Codes.md"
    	],
		"Examples" => [
			"Examples/The Vardy-Be’ery Decomposition.md",
			"Examples/Quantum Reed-Muller Codes.md"
		],
		"Classical" => [
        	"Classical/linear_code.md",
			"Classical/concatenation.md",
        	"Classical/cyclic_code.md",
        	"Classical/quasi-cyclic_code.md",
        	"Classical/GeneralizedReedSolomon.md",
        	"Classical/ReedMuller.md",
			"Classical/new_codes_from_old.md",
			"Classical/product_codes.md",
        	"Classical/misc_known_codes.md"
    	],
		"LDPC" => [
			"LDPC/codes.md",
			"LDPC/Tanner_codes.md",
			"LDPC/analysis.md",
			"LDPC/channels.md",
			"LDPC/decoders.md"
		],
    	"Quantum" => [
        	"Quantum/quantum_code.md",
			"Quantum/product_codes.md",
        	"Quantum/misc_known_codes.md"
    	],
    	"Misc" => [
        	"tilings.md",
        	"trellis.md",
        	"utils.md",
        	"weight_dist.md"
    	],
    	"References" => "references.md",
    	"Index" => "theindex.md"
    # "Developer Documentation" => [

    # ],
	]
)

deploydocs(repo = "github.com/esabo/CodingTheory.git", devbranch = "subsystem")
