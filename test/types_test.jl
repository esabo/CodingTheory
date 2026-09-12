@testitem "Code type hierarchy" begin
    using Oscar, CodingTheory

    @testset "Classical Supertypes" begin
        for T in (AbstractNonadditiveCode, AbstractNonlinearCode, AbstractAdditiveCode)
            @test T <: AbstractCode
        end
        @test AbstractLinearCode <: AbstractAdditiveCode

        for T in (
            AbstractMatrixProductCode,
            AbstractReedMullerCode,
            AbstractCyclicCode,
            AbstractQuasiCyclicCode,
            AbstractGeneralizedReedSolomonCode,
            AbstractAlgebraicGeometryCode,
            AbstractConcatenatedCode,
            AbstractAlternateCode,
            AbstractTwistedReedSolomonCode,
            AbstractTannerCode,
        )
            @test T <: AbstractLinearCode
        end

        # cyclic codes specialize to BCH and then to Reed--Solomon
        @test AbstractBCHCode <: AbstractCyclicCode
        @test AbstractReedSolomonCode <: AbstractBCHCode
        @test AbstractCyclicCode2D <: AbstractCyclicCode

        # Goppa and generalized Srivastava codes are alternate codes
        @test AbstractGoppaCode <: AbstractAlternateCode
        @test AbstractGeneralizedSrivastavaCode <: AbstractAlternateCode

        # the hierarchy is not accidentally collapsed
        @test !(AbstractLinearCode <: AbstractNonlinearCode)
        @test !(AbstractCyclicCode <: AbstractBCHCode)
    end

    @testset "Concrete Classical Codes" begin
        @test HammingCode(2, 3) isa AbstractLinearCode
        @test BCHCode(2, 15, 3) isa AbstractBCHCode
        @test ReedSolomonCode(8, 3) isa AbstractReedSolomonCode
        @test ReedMullerCode(1, 3) isa AbstractLinearCode
        @test LDPCCode(parity_check_matrix(HammingCode(2, 3))) isa AbstractLDPCCode
    end

    @testset "LDPC Channels" begin
        @test AbstractLDPCCode <: AbstractLinearCode
        # AbstractChannel must be qualified because Base exports the same name
        @test AbstractDiscreteChannel <: CodingTheory.AbstractChannel
        @test AbstractContinuousChannel <: CodingTheory.AbstractChannel
        @test BinarySymmetricChannel(0.1) isa AbstractDiscreteChannel
        @test BinaryErasureChannel(0.1) isa AbstractDiscreteChannel
        @test BAWGNChannel(2.0) isa AbstractContinuousChannel

        # the exported abbreviations name the same types
        @test BSC === BinarySymmetricChannel
        @test BEC === BinaryErasureChannel
        @test BAWGNC === BAWGNChannel
    end

    @testset "Quantum Supertypes" begin
        @test AbstractSubsystemCode <: AbstractAdditiveCode
        @test AbstractStabilizerCode <: AbstractSubsystemCode
        @test AbstractSubsystemCodeCSS <: AbstractSubsystemCode
        @test AbstractStabilizerCodeCSS <: AbstractStabilizerCode
        @test AbstractHypergraphProductCode <: AbstractStabilizerCodeCSS

        for T in (AbstractGraphStateSubsystem, AbstractGraphStateStabilizer)
            @test T <: AbstractSubsystemCode
        end
        @test AbstractGraphStateSubsystemCSS <: AbstractSubsystemCodeCSS
        @test AbstractGraphStateStabilizerCSS <: AbstractStabilizerCodeCSS

        # entanglement-assisted codes mirror the unassisted hierarchy
        @test AbstractEASubsystemCode <: AbstractSubsystemCode
        @test AbstractEASubsystemCodeCSS <: AbstractEASubsystemCode
        @test AbstractEAStabilizerCode <: AbstractStabilizerCode
        @test AbstractEAStabilizerCodeCSS <: AbstractEAStabilizerCode

        # a stabilizer code is a subsystem code, but not conversely
        @test !(AbstractSubsystemCode <: AbstractStabilizerCode)
    end
end

@testitem "Holy traits" begin
    using Oscar, CodingTheory

    # dispatch throughout the quantum code is driven by these three traits
    @testset "Trait Types" begin
        @test HasLogicals <: LogicalTrait
        @test HasNoLogicals <: LogicalTrait
        @test HasGauges <: GaugeTrait
        @test HasNoGauges <: GaugeTrait
        @test IsCSS <: CSSTrait
        @test IsNotCSS <: CSSTrait
    end

    @testset "Stabilizer Codes" begin
        S = SteaneCode()
        T = typeof(S)

        @test LogicalTrait(T) == HasLogicals()
        # a stabilizer code has no gauge degrees of freedom
        @test GaugeTrait(T) == HasNoGauges()
        @test CSSTrait(T) == IsCSS()
        @test is_CSS(S)

        # the trait agrees with the CSS-ness of the code
        NC = Q513()
        @test CSSTrait(typeof(NC)) == IsNotCSS()
        @test !is_CSS(NC)
        @test GaugeTrait(typeof(NC)) == HasNoGauges()
    end

    @testset "Subsystem Codes" begin
        S = GaugedShorCode()
        T = typeof(S)

        @test GaugeTrait(T) == HasGauges()
        @test LogicalTrait(T) == HasLogicals()
        # gauge operators exist beyond the stabilizers
        @test rank(gauge_group(S)) > rank(stabilizers(S))
    end

    @testset "Graph States" begin
        # the graph-state types are the only ones carrying no logicals
        for T in (
            CodingTheory.GraphStateStabilizer,
            CodingTheory.GraphStateStabilizerCSS,
            CodingTheory.GraphStateSubsystem,
            CodingTheory.GraphStateSubsystemCSS,
        )
            @test LogicalTrait(T) == HasNoLogicals()
        end
        @test GaugeTrait(CodingTheory.GraphStateSubsystem) == HasGauges()
        @test GaugeTrait(CodingTheory.GraphStateStabilizer) == HasNoGauges()

        # a stabilizer group of full rank encodes nothing, though the ordinary
        # constructor still returns a StabilizerCode
        F = Oscar.Nemo.Native.GF(2)
        S = StabilizerCode(matrix(F, [1 0 0 1; 0 1 1 0]))
        @test S.k == 0
    end
end
