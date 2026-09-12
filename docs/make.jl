using DocumenterCitations, Documenter, CodingTheory

bib = CitationBibliography(joinpath(@__DIR__, "src", "references.bib"); style = :numeric)

Documenter.makedocs(;
    plugins = [bib],
    clean = true,
    doctest = false,
    modules = Module[CodingTheory],
    repo = "https://github.com/esabo/CodingTheory/blob/{commit}{path}#L{line}",
    highlightsig = true,
    sitename = "CodingTheory.jl",
    warnonly = [:missing_docs],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        collapselevel = 1,
    ),
    pages = [
        "Home" => "index.md",
        "Getting started" => [
            "Linear codes" => "Tutorials/Linear Codes.md",
            "Cyclic codes" => "Tutorials/Cyclic Codes.md",
            "Quantum codes" => "Tutorials/Quantum Codes.md",
            "Message-passing decoding" => "Tutorials/Message Passing.md",
            "Minimum distance" => "Tutorials/min_dist.md",
            "Input and output" => "Tutorials/Input and Output.md",
            "Weight reduction" => "Tutorials/Weight Reduction.md",
        ],
        "Classical codes" => [
            "Core API" => "Classical/linear_code.md",
            "Cyclic codes" => "Classical/cyclic_code.md",
            "Quasi-cyclic codes" => "Classical/quasi-cyclic_code.md",
            "Generalized Reed--Solomon codes" => "Classical/GeneralizedReedSolomon.md",
            "Reed--Muller codes" => "Classical/ReedMuller.md",
            "New codes from old" => "Classical/new_codes_from_old.md",
            "Concatenation" => "Classical/concatenation.md",
            "Product codes" => "Classical/product_codes.md",
            "Known codes" => "Classical/misc_known_codes.md",
            "Minimum distance" => "Classical/min_dist.md",
            "Information-set decoding" => "Classical/ISD_attacks.md",
            "Designs and invariant theory" => "Classical/designs.md",
        ],
        "LDPC codes" => [
            "Constructions" => "LDPC/codes.md",
            "Tanner codes" => "LDPC/Tanner_codes.md",
            "Analysis" => "LDPC/analysis.md",
            "Cycles and ACE" => "LDPC/cycles.md",
            "Channels" => "LDPC/channels.md",
            "Decoders" => "LDPC/decoders.md",
            "Generalized belief propagation" => "LDPC/GBP.md",
        ],
        "Quantum codes" => [
            "Core API" => "Quantum/quantum_code.md",
            "Subsystem codes" => "Quantum/subsystem_code.md",
            "Product and BB codes" => "Quantum/product_codes.md",
            "Known codes" => "Quantum/misc_known_codes.md",
            "New codes from old" => "Quantum/new_codes_from_old.md",
            "Minimum distance" => "Quantum/min_dist.md",
            "Bounds and weight enumerators" => "Quantum/bounds.md",
            "Quantum LDPC parameters" => "Quantum/QLDPC.md",
            "Weight reduction" => "Quantum/weight_reduction.md",
            "Homological measurements" => "Quantum/homological_measurements.md",
            "Expansion and confinement" => "Quantum/code_expansion.md",
        ],
        "Algorithms and utilities" => [
            "Weight enumerators" => "weight_dist.md",
            "Trellises" => "trellis.md",
            "Tilings" => "tilings.md",
            "Input and output" => "io.md",
            "Utilities" => "utils.md",
        ],
        "Examples" => [
            "Vardy--Be'ery decomposition" => "Examples/The Vardy-Be’ery Decomposition.md",
            "Quantum Reed--Muller codes" => "Examples/Quantum Reed-Muller Codes.md",
        ],
        "References" => "references.md",
        "Function index" => "theindex.md",
    ],
)

deploydocs(repo = "github.com/esabo/CodingTheory.git", devbranch = "dev")
