using UncertainHistogramming
using Test

@testset "UncertainHistogramming.jl" begin
    
    @time @testset "Single values" begin
        μ, σ = 1.0, 4.0
        ghist = GaussianHistogram()
        uhist = UniformHistogram()
        push!(ghist, (μ, σ))
        push!(uhist, (μ, σ))

        @time @testset "1-Point Structure" begin
            @test length(ghist) == length(ghist.values)
            @test length(ghist) == length(ghist.errors)
            @test length(uhist) == length(uhist.values)
            @test length(uhist) == length(uhist.errors)
        end

        @time @testset "Gaussian 1-Point Statistics" begin 
            @test mean(ghist) == μ
            @test var(ghist) == σ^2
            @test std(ghist) == σ 
        end
        
        @time @testset "Uniform 1-Point Statistics" begin 
            @test mean(uhist) == μ
            @test var(uhist) ≈ σ^2 / 3
            @test std(uhist) ≈ σ / sqrt(3)
        end

        μ1, σ1 = -1.0, 0.5
        push!(ghist, (μ1, σ1))
        push!(uhist, (μ1, σ1))

        @time @testset "2-Point Structure" begin
            @test length(ghist) == length(ghist.values)
            @test length(ghist) == length(ghist.errors)
            @test length(uhist) == length(uhist.values)
            @test length(uhist) == length(uhist.errors)
        end

        @time @testset "Gaussian 2-Point Statistics" begin 
            @test mean(ghist) == 0.5 * ( μ + μ1 )
            @test var(ghist) == 0.5 * ( σ^2 + σ1^2 ) + 0.25 * ( μ - μ1 )^2
            @test std(ghist) == sqrt( 0.5 * ( σ^2 + σ1^2 ) + 0.25 * ( μ - μ1 )^2 ) 
        end
        
        @time @testset "Uniform 2-Point Statistics" begin 
            @test mean(uhist) == 0.5 * ( μ + μ1 )
            @test var(uhist) ≈ 0.5 * ( (σ^2 + σ1^2) / 3 ) + 0.25 * ( μ - μ1 )^2
            @test std(uhist) ≈ sqrt( 0.5 * ( (σ^2 + σ1^2) / 3 ) + 0.25 * ( μ - μ1 )^2 ) 
        end
    end

end
