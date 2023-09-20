# TODO: fix
# @testset "QWEMacWId" begin
#     using CodingTheory

#     # basic - too symmetric in X, Y, Z to determine errors in formula
#     Q = SteaneCode()
#     WEstabs = CodingTheory._weightenumeratorBFQ(Q.stabs, Q.charvec, missing)
#     WEnorm = CodingTheory._weightenumeratorBFQ(Q.dualgens, Q.charvec, parent(WEstabs.polynomial))
#     WEnormMacW = MacWilliams_identity(Q, WEstabs)
#     WEstabsMacW = MacWilliams_identity(Q, WEnorm, true)
#     @test WEnormMacW == WEnorm
#     @test WEstabsMacW == WEstabs

#     # non-CSS - also too symmetric in X, Y, Z
#     Q = Q513()
#     WEstabs = CodingTheory._weightenumeratorBFQ(Q.stabs, Q.charvec, missing)
#     WEnorm = CodingTheory._weightenumeratorBFQ(Q.dualgens, Q.charvec, parent(WEstabs.polynomial))
#     WEnormMacW = MacWilliams_identity(Q, WEstabs)
#     WEstabsMacW = MacWilliams_identity(Q, WEnorm, true)
#     @test WEnormMacW == WEnorm
#     @test WEstabsMacW == WEstabs

#     # k > 1 - magically also too symmetric in X, Y, Z
#     Q = Q823()
#     WEstabs = CodingTheory._weightenumeratorBFQ(Q.stabs, Q.charvec, missing)
#     WEnorm = CodingTheory._weightenumeratorBFQ(Q.dualgens, Q.charvec, parent(WEstabs.polynomial))
#     WEnormMacW = MacWilliams_identity(Q, WEstabs)
#     WEstabsMacW = MacWilliams_identity(Q, WEnorm, true)
#     @test WEnormMacW == WEnorm
#     @test WEstabsMacW == WEstabs

#     # # k > 1
#     # Q = Q1573()
#     # WEstabs = CodingTheory._weightenumeratorBFQ(Q.stabs, Q.charvec, missing)
#     # WEnorm = CodingTheory._weightenumeratorBFQ(Q.dualgens, Q.charvec, parent(WEstabs.polynomial))
#     # WEnormMacW = MacWilliams_identity(Q, WEstabs)
#     # WEstabsMacW = MacWilliams_identity(Q, WEnorm, true)
#     # @test WEnormMacW == WEnorm
#     # @test WEstabsMacW == WEstabs

#     # popular - non-symmetric, detected error in X and Z terms being switched in MacWilliams_identity
#     Q = Q15RM()
#     WEstabs = CodingTheory._weightenumeratorBFQ(Q.stabs, Q.charvec, missing)
#     WEnorm = CodingTheory._weightenumeratorBFQ(Q.dualgens, Q.charvec, parent(WEstabs.polynomial))
#     WEnormMacW = MacWilliams_identity(Q, WEstabs)
#     WEstabsMacW = MacWilliams_identity(Q, WEnorm, true)
#     @test WEnormMacW == WEnorm
#     @test WEstabsMacW == WEstabs
# end
