# Reproduce the LP₂ code from Table I of arXiv:2410.02753.
using CodingTheory
using Dates
using GLPK
using JuMP
using Oscar

function lp2_code()
    F = Oscar.Nemo.Native.GF(2)
    S, x = polynomial_ring(F, :x)
    R, _ = residue_ring(S, x^9 - 1)
    A = matrix(R, 3, 4, [
        1, 1, 1, 1,
        1, x, x^6, x^7,
        1, x^4, x^5, x^2,
    ])
    B = CodingTheory._CT_adjoint(A)

    I3 = identity_matrix(R, 3)
    I4 = identity_matrix(R, 4)
    H_X = lift(hcat(
        kronecker_product(A, I4),
        kronecker_product(I3, B),
    ))
    H_Z = lift(hcat(
        kronecker_product(I4, CodingTheory._CT_adjoint(B)),
        kronecker_product(CodingTheory._CT_adjoint(A), I3),
    ))
    return CSSCode(H_X, H_Z)
end

max_d = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 12
time_limit_sec = if length(ARGS) < 2 || lowercase(ARGS[2]) == "none"
    nothing
else
    parse(Float64, ARGS[2])
end

println("[$(now())] Constructing LP₂ parity-check matrices")
flush(stdout)
code = lp2_code()

println("parameters = [[$(length(code)), $(dimension(code))]]")
println("[$(now())] Computing the logical basis")
flush(stdout)
logicals_matrix(code)
println("[$(now())] Logical basis ready")
flush(stdout)
for which in (:X, :Z)
    println("[$(now())] Starting exact $which solve through weight $max_d")
    flush(stdout)
    result = nothing
    elapsed = @elapsed result = minimum_distance(
        code; which=which, alg=:ILP, max_d=max_d,
        time_limit_sec=time_limit_sec, verbose=true,
    )
    d, witness = result
    println("[$(now())] Finished exact $which solve")
    println("$which distance through $max_d = $d")
    println("$which elapsed seconds = $elapsed")
    if d != -1
        println("$which witness support = $(findall(!iszero, collect(witness)))")
        println("$which witness logical = $(is_logical(code, witness))")
    end
    flush(stdout)
end
