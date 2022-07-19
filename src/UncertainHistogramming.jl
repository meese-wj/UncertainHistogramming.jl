module UncertainHistogramming

import Base: eltype, push!, length, getindex, setindex!, show
import StaticArrays: MArray, @MArray, @SArray
import Statistics: mean, var, std
import StatsBase: skewness, kurtosis
import Measurements: Measurement, measurement
using RecipesBase

include("Moments.jl")

export 
# Base overloads in ./util.jl
       push!, eltype, length, getindex, setindex!, show,
# Statistics and StatsBase overloads in ./stats.jl
       mean, var, std, skewness, kurtosis,
# Measurements.jl overloads
       measurement,
# ./Moments.jl overloads
       moment, FirstMoment, SecondMoment, ThirdMoment, FourthMoment,
# UncertainHistogramming exports
       ContinuousHistogram, GaussianHistogram, 
       gaussian, construct, construct!, kernel, val_err

"""
    abstract type ContinuousHistogram end

`abstract type` representing the `ContinuousHistogram` concept. These are
histograms whose input data have _uncertainty_ associated with them, and therefore
we build them using a value-error-dependent kernel for each entry.
"""
abstract type ContinuousHistogram end
@doc raw"""
    kernel(::ContinuousHistogram, ::AbstractArray, data) -> MethodError
    kernel(::GaussianHistogram, ::AbstractArray, data)

The kernel function ``\mathcal{K}(x)`` used to compute a [`ContinuousHistogram`](@ref) such that

```math
\mathcal{H}(x) = \sum_{i = 1}^M \mathcal{K}_i(x),
```

where ``M`` is the total number of data points. For brevity, we have represented all `data`-dependence
through the subscript ``i``. 

!!! note
    By default, we `MethodError` out for any `<: ContinuousHistogram` until its directly implemented.
"""
kernel(hist::ContinuousHistogram, args...) = throw( MethodError(kernel, args...) )
construct(hist::ContinuousHistogram, x) = throw( MethodError(construct, hist, x) )
construct!(output, hist::ContinuousHistogram, x) = throw( MethodError(construct!, output, hist, x) )


@doc raw"""
    GaussianHistogram{T <: Number} <: ContinuousHistogram

A [`ContinuousHistogram`](@ref) with a Gaussian kernel for each value-error pair.

This [`ContinuousHistogram`] is formed by summing Gaussian kernels for each ``(\mu_i, \sigma_i)`` as

```math
\mathcal{H}(y) = \frac{1}{M} \sum_{i = 1}^M \frac{ \exp\left[ -\frac{ \left( y - \mu_i \right)^2 }{2 \sigma_i^2} \right] }{ \sigma_i \sqrt{2\pi}}.
```

This expression makes the calculation of the non-central moments a simple mean of the 
individual non-central moments.

# Contents
- `moments::Vector{T}`: a collection of moments for the `GaussianHistogram` which are updated in an _online_ fashion
- `values::Vector{T}`: the values used to [`construct`](@ref) the `GaussianHistogram`
- `errors::Vector{T}`: the errors used to [`construct`](@ref) the `GaussianHistogram`

!!! note
    The statistics come from calculations involving the `moments`. The `values` and 
    `errors` are necessarily stored for visualization purposes.
"""
mutable struct GaussianHistogram{T <: Number} <: ContinuousHistogram
    moments::Vector{T}
    values::Vector{T}
    errors::Vector{T}

    GaussianHistogram{T}() where {T} = ( temp = @MArray zeros(T, 4); new( temp, zeros(T, 0), zeros(T, 0) ) )
    GaussianHistogram(args...) = GaussianHistogram{Float64}(args...)
end

include("util.jl")
include("stats.jl")
include("plotrecipes.jl")

"""
    gaussian(::Number, μ, σ)
    gaussian(::AbstractArray, μ, σ)

Calculate the normalized value of a Gaussian with mean μ and variance σ².
"""
gaussian(x::Number, μ, σ) = exp( -0.5 * (x - μ)^2 / σ^2 ) / ( σ * sqrt(2π) )
gaussian(x::AbstractArray, μ, σ) = broadcast( y -> gaussian(y, μ, σ), x )
kernel(hist::GaussianHistogram, x::AbstractArray, data) where {T} = gaussian(x, data...) / length(hist)

"""
    val_err(::GaussianHistogram, idx)

Return a `(value, error)`-`Tuple`.
"""
val_err(hist::GaussianHistogram, idx) = (hist.values[idx], hist.errors[idx])

"""
    construct!(output, ::GaussianHistogram, x)

Similar to [`construct`](@ref) but here, the `output` `Array` is modified in-place.
"""
function construct!(output, hist::GaussianHistogram, x)
    # TODO: I think I can define this just for a ContinuousHistogram only as long as 
    #       I demand each ContinuousHistogram has a val_err method.
    size(x) == size(output) ? nothing : ArgumentError("input and output vectors are of different sizes: $(size(x)) != $(size(output))")
    for (μ, σ) ∈ zip(hist.values, hist.errors)
        @views output .+= kernel(hist, x, (μ, σ))
    end
    return output
end

"""
    construct(::GaussianHistogram, x)

Map the values of `x` through the [`GaussianHistogram`] and return an
`Array` of the same `size` as `x`.
"""
function construct(hist::GaussianHistogram, x)
    output = zeros(size(x))
    construct!(output, hist, x)
end

"""
    measurement(::GaussianHistogram)

Interface to [`Measurements.jl`](https://juliaphysics.github.io/Measurements.jl/stable/). 
Return a `Measurement` with `val` as the [`GaussianHistogram`](@ref) [`mean`](@ref) and the 
`err` as the [`GaussianHistogram`](@ref) [`std`](@ref) (standard deviation).
"""
measurement(hist::GaussianHistogram) = measurement(mean(hist), std(hist))

end
