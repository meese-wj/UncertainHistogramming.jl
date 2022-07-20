
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

"""
    gaussian(::Number, μ, σ)
    gaussian(::AbstractArray, μ, σ)

Calculate the normalized value of a Gaussian with mean μ and variance σ².
"""
gaussian(x::Number, μ, σ) = exp( -0.5 * (x - μ)^2 / σ^2 ) / ( σ * sqrt(2π) )
gaussian(x::AbstractArray, μ, σ) = broadcast( y -> gaussian(y, μ, σ), x )
kernel(hist::GaussianHistogram, x::AbstractArray, data) where {T} = gaussian(x, data...) / length(hist)
