using UncertainHistogramming
using Test

@testset "UncertainHistogramming.jl" begin
    
    @time @testset "Single values" begin
        μ, σ = 1.0, 4.0
        hist = GaussianHistogram()
        push!(hist, (μ, σ))

        @time @testset "1-Point Structure" begin
            @test length(hist) == length(hist.values)
            @test length(hist) == length(hist.errors)
        end

        @time @testset "1-Point Statistics" begin 
            @test mean(hist) == μ
            @test var(hist) == σ^2
            @test std(hist) == σ 
        end

        μ1, σ1 = -1.0, 0.5
        push!(hist, (μ1, σ1))

        @time @testset "2-Point Structure" begin
            @test length(hist) == length(hist.values)
            @test length(hist) == length(hist.errors)
        end

        @time @testset "2-Point Statistics" begin 
            @test mean(hist) == 0.5 * ( μ + μ1 )
            @test var(hist) == 0.5 * ( σ^2 + σ1^2 ) + 0.25 * ( μ - μ1 )^2
            @test std(hist) == sqrt( 0.5 * ( σ^2 + σ1^2 ) + 0.25 * ( μ - μ1 )^2 ) 
        end
    end

end
