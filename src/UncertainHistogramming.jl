module UncertainHistogramming

import Base: eltype, push!, length, getindex, setindex!, show
import StaticArrays: MArray, @MArray
import Statistics: mean, var, std
import Measurements: Measurement, measurement

include("Moments.jl")

export 
# Base overloads in ./util.jl
       push!, eltype, length, getindex, setindex!, show,
# Statistics overloads in ./stats.jl
       mean, var, std,
# Measurements.jl overloads
       measurement,
# ./Moments.jl overloads
       moment, skewness, kurtosis,
# UncertainHistogramming exports
       ContinuousHistogram, GaussianHistogram, gaussian

abstract type ContinuousHistogram end
construct(hist::ContinuousHistogram, x) = throw( MethodError(construct, hist, x) )
construct!(output, hist::ContinuousHistogram, x) = throw( MethodError(construct!, output, hist, x) )

mutable struct GaussianHistogram{T <: Number}
    moments::Vector{T}
    values::Vector{T}
    errors::Vector{T}

    GaussianHistogram{T}() where {T} = new( ( zeros(T, 4) ), zeros(T, 0), zeros(T, 0) )
    GaussianHistogram(args...) = GaussianHistogram{Float64}(args...)
end

include("util.jl")
include("stats.jl")

gaussian(x::Number, μ, σ) = exp( -0.5 * (x - μ)^2 / σ^2 ) / ( σ * sqrt(2π) )
gaussian(x::AbstractArray, μ, σ) = broadcast( y -> gaussian(y, μ, σ), x )

function construct!(output, hist::GaussianHistogram, x)
    size(x) == size(output) ? nothing : ArgumentError("input and output vectors are of different sizes: $(size(x)) != $(size(output))")
    for (μ, σ) ∈ zip(hist.values, hist.errors)
        @views output .+= gaussian(x, μ, σ)
    end
    return output ./ length(hist)
end

function construct(hist::GaussianHistogram, x)
    output = zeros(size(x))
    construct!(output, hist, x)
end

measurement(hist::GaussianHistogram) = measurement(mean(hist), std(hist))

end
