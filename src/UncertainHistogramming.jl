module UncertainHistogramming

import Base: eltype, push!, length, getindex, setindex!, show
import StaticArrays: MArray, @MArray
import Statistics: mean, var, std
import Measurements: Measurement, measurement
include("Moments.jl")

export 
# Base overloads
       push!, eltype, length, getindex, setindex!, show,
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
    moments::Vector{T}
    values::Vector{T}
    errors::Vector{T}

    GaussianHistogram{T}() where {T} = new( ( zeros(T, 4) ), zeros(T, 0), zeros(T, 0) )
    GaussianHistogram(args...) = GaussianHistogram{Float64}(args...)
end

gaussian(x::Number, μ, σ) = exp( -0.5 * (x - μ)^2 / σ^2 ) / ( σ * sqrt(2π) )
gaussian(x::AbstractArray, μ, σ) = broadcast( y -> gaussian(y, μ, σ), x )

eltype(::GaussianHistogram{T}) where {T} = T
length(hist::GaussianHistogram) = length(hist.values)

mean( hist::GaussianHistogram ) = moment(hist, FirstMoment)
var( hist::GaussianHistogram ) = moment(hist, SecondMoment) - mean(hist)^2
std( hist::GaussianHistogram ) = sqrt( var(hist) )

function push!( hist::GaussianHistogram, tup::Tuple{Number, Number} )
    _update_moments!(hist, tup...)
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

for moment_t ∈ moment_list
    @eval Base.getindex(hist::GaussianHistogram, ::Type{$moment_t}) = hist.moments[ MomentIndex($moment_t) ]
    @eval Base.setindex!(hist::GaussianHistogram{T}, val::S, ::Type{$moment_t}) where {T, S} = hist.moments[ MomentIndex($moment_t) ] = convert(T, val)
    @eval moment(hist::GaussianHistogram, ::Type{$moment_t}) = hist[$moment_t]
end

measurement(hist::GaussianHistogram) = measurement(mean(hist), std(hist))

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

"""
    _update_moments!(::GaussianHistogram, μ, σ)

Update the the non-central moments of the [`GaussianHistogram`](@ref) _online_ with
the new value-error pair `(μ, σ)`. This is done because each is non-central moment
is the mean non-central moments of each [`GaussianDistribution`](@ref).
"""
function _update_moments!(hist::GaussianHistogram, μ, σ)
    _update_moment!(hist, FirstMoment, μ, σ)
    _update_moment!(hist, SecondMoment, μ, σ)
    _update_moment!(hist, ThirdMoment, μ, σ)
    _update_moment!(hist, FourthMoment, μ, σ)
    return nothing
end

"""
    _update_moment!(::GaussianHistogram, moment_t, μ, σ)

Update the [`GaussianHistogram`](@ref)'s `moment_t` using an [`_online_mean`](@ref)
with the inclusion of the value-error pair `(μ, σ)`.
"""
@inline function _update_moment!(hist::GaussianHistogram, moment_t, μ, σ) 
    hist[moment_t] = _online_mean( moment(GaussianDistribution, moment_t, μ, σ), moment(hist, moment_t), length(hist) )
    return nothing
end

@doc raw"""
    _online_mean(xnew, μold, nold)

Update the old mean `μold` over `nold` elements with the new data point `xnew`.
This function implements

```math
\mu_{n+1} = \mu_n + \frac{x_{n+1} - \mu_{n}}{n+1}.
```
"""
@inline _online_mean(xnew, μold, nold) = μold + (xnew - μold) / (nold + 1)

"""
    show([::IO,] ::GaussianHistogram)
    show(::GaussianHistogram)

`print` the relevant information for a [`GaussianHistogram`](@ref).
"""
function show(io::IO, hist::GaussianHistogram)
    println(io, "GaussianHistogram{$(eltype(hist))}:")
    println(io, "  length  = $(length(hist))")
    momstring = ""
    for mom ∈ hist.moments momstring *= "$mom  " end
    println(io, "  moments = $(momstring)")
    println(io, "")
    println(io, "  Statistics")
    println(io, "    mean        = $(mean(hist))")
    println(io, "    variance    = $(var(hist))")
    println(io, "    skewness    = $(skewness(hist))")
    println(io, "    kurtosis    = $(kurtosis(hist))")
end
show(hist::GaussianHistogram) = show(stdout, hist)

end
