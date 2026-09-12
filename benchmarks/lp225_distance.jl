# Reproduce the LP₂ code from Table I of arXiv:2410.02753.
# Default: X-only HiGHS benchmark through weight 12, with the order-9
# lifted circulant symmetry cut and 4 solver threads (the performance
# cores on this 4P+6E machine).
using CodingTheory
using Dates
using HiGHS
using JuMP

const Oscar = CodingTheory.Oscar

function lp2_code()
    F = Oscar.Nemo.Native.GF(2)
    S, x = Oscar.polynomial_ring(F, :x)
    R, _ = Oscar.residue_ring(S, x^9 - 1)
    A = Oscar.matrix(R, 3, 4, [
        1, 1, 1, 1,
        1, x, x^6, x^7,
        1, x^4, x^5, x^2,
    ])
    B = CodingTheory._CT_adjoint(A)

    I3 = Oscar.identity_matrix(R, 3)
    I4 = Oscar.identity_matrix(R, 4)
    H_X = CodingTheory.lift(hcat(
        kronecker_product(A, I4),
        kronecker_product(I3, B),
    ))
    H_Z = CodingTheory.lift(hcat(
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
ilp_optimizer = length(ARGS) >= 3 ? Symbol(ARGS[3]) : :HiGHS
ilp_threads = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 4
ilp_cyclic_period = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 9

println("[$(now())] Constructing LP₂ parity-check matrices")
flush(stdout)
code = lp2_code()

println("parameters = [[$(length(code)), $(dimension(code))]]")
println("[$(now())] Computing the logical basis")
flush(stdout)
logicals_matrix(code)
println("[$(now())] Logical basis ready")
println("[$(now())] HiGHS X-only benchmark: max_d=$max_d, threads=$ilp_threads, cyclic_period=$ilp_cyclic_period")
flush(stdout)

println("[$(now())] Starting exact X solve through weight $max_d")
flush(stdout)
result = nothing
elapsed = @elapsed result = minimum_distance(
    code; which=:X, alg=:ILP, max_d=max_d,
    time_limit_sec=time_limit_sec, ilp_optimizer=ilp_optimizer,
    ilp_threads=ilp_threads, ilp_cyclic_period=ilp_cyclic_period,
    verbose=true,
)
d, witness = result
println("[$(now())] Finished exact X solve")
println("X distance through $max_d = $d")
println("X elapsed seconds = $elapsed")
if d != -1
    println("X witness support = $(findall(!iszero, collect(witness)))")
    println("X witness logical = $(is_logical(code, witness))")
end
flush(stdout)
