module UncertainHistogramming

import Base: eltype, push!, length
import Statistics: mean, var, std
import Measurements: Measurement, measurement
import Moments: GaussianDistribution,
                moment, FirstMoment, SecondMoment, ThirdMoment, FourthMoment
                skewness, kurtosis

export 
# Base overloads
       push!, eltype, length,
# Statistics overloads
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
    values::Vector{T}
    errors::Vector{T}

    GaussianHistogram{T}() where {T} = new( zeros(T, 0), zeros(T, 0) )
    GaussianHistogram(args...) = GaussianHistogram{Float64}(args...)
end

gaussian(x::Number, μ, σ) = exp( -0.5 * (x - μ)^2 / σ^2 ) / ( σ * sqrt(2π) )
gaussian(x::AbstractArray, μ, σ) = broadcast( y -> gaussian(y, μ, σ), x )

eltype(::GaussianHistogram{T}) where {T} = T
length(hist::GaussianHistogram) = length(hist.values)

mean( hist::GaussianHistogram ) = mean(hist.values)
var( hist::GaussianHistogram ) = mean(x -> x^2, hist.errors) + var(hist.values, corrected = false)
std( hist::GaussianHistogram ) = sqrt( var(hist) )

function push!( hist::GaussianHistogram, tup::Tuple{Number, Number} )
    push!(hist.values, tup[1])
    push!(hist.errors, tup[2])
    return hist
end
push!(hist::GaussianHistogram, meas::Measurement) = push!(hist, (meas.val, meas.err))

function push!(hist::GaussianHistogram, tup_vec::Vector{Tuple{Number, Number}})
    for tup ∈ tup_vec
        push!(hist, tup)
    end
    return hist
end
function push!(hist::GaussianHistogram, meas_vec::Vector{Measurement})
    for meas ∈ meas_vec
        push!(hist, meas)
    end
    return hist
end

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

moment(hist::GaussianHistogram, ::Type{FirstMoment}) = mean(hist)
function _moment(hist::GaussianHistogram, moment_t)
    val = zero(eltype(hist))
    for (μ, σ) ∈ zip(hist.values, hist.errors)
        val += moment(GaussianDistribution, moment_t, μ, σ)
    end
    return val / length(hist)
end

for moment_t ∈ (:SecondMoment, :ThirdMoment, :FourthMoment)
    @eval moment(hist::GaussianHistogram, ::Type{$moment_t}) = _moment(hist, $moment_t)
end

"""
    skewness(::GaussianHistogram)

[Fisher's skewness](https://en.wikipedia.org/wiki/Skewness?oldformat=true#Fisher's_moment_coefficient_of_skewness) of a [`GaussianHistogram`](@ref).

!!! note
    For comparison purposes, this `skewness` definition applied to a Gaussian distribution yields identically zero.
"""
function skewness(hist::GaussianHistogram)
    μ = mean(hist)
    val =  moment(hist, ThirdMoment) - 3 * μ * moment(hist, SecondMoment) + 2 * μ^3
    return val / std(hist)^3
end

"""
    kurtosis(::GaussianHistogram [, excess = false ])

[Pearson kurtosis](https://en.wikipedia.org/wiki/Kurtosis?oldformat=true#Pearson_moments) of a [`GaussianHistogram`](@ref).

!!! note
    For comparison purposes, this `kurtosis` definition, with `excess == false`, applied to a 
    Gaussian distribution yields `3`. If `excess == true`, then the `kurtosis` vanishes for a Gaussian.
"""
function kurtosis(hist::GaussianHistogram, excess = false)
    μ = mean(hist)
    val =  moment(hist, FourthMoment) - 4 * μ * moment(hist, ThirdMoment)
    val += 6 * μ^2 * moment(hist, SecondMoment) - 3 * μ^4
    val /= var(hist)^2
    return excess ? val - convert(typeof(val), 3) : val
end



end
