module UncertainHistogramming

import Base: eltype, push!, length
import Statistics: mean, var, std
import Measurements: Measurement, measurement

export 
# Base overloads
       push!, eltype, length,
# Statistics overloads
       mean, var, std,
# Measurements.jl overloads
       measurement,
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

end
