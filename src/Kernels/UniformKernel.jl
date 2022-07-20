
@doc raw"""
    UniformHistogram{T <: Number} <: ContinuousHistogram

A [`ContinuousHistogram`](@ref) with a uniform kernel for each value-error pair.

This [`ContinuousHistogram`] is formed by summing uniform kernels for each ``(x_i, \epsilon_i)`` as

```math
\mathcal{H}(y) = \frac{1}{M} \sum_{i = 1}^M \mathcal{U}(y; x_i, \epsilon_i),
```

where 

```math
\mathcal{U}(y; x_i, \epsilon_i) = \begin{cases}
\frac{1}{2\epsilon_i}, & y \in (x_i - \epsilon_i, x_i + \epsilon_i)
\\
0, & \mathrm{otherwise}
\end{cases}.
```

These expression makes the calculation of the non-central [`moment`](@ref)s a simple mean of the 
individual non-central [`moment`](@ref).

# Contents
- `moments::Vector{T}`: a collection of moments for the `UniformHistogram` which are updated in an _online_ fashion
- `values::Vector{T}`: the values used to [`construct`](@ref) the `UniformHistogram`
- `errors::Vector{T}`: the errors used to [`construct`](@ref) the `UniformHistogram`

!!! note
    The statistics come from calculations involving the [`moment`](@ref). The `values` and 
    `errors` are necessarily stored for visualization purposes.
"""
mutable struct UniformHistogram{T <: Number} <: ContinuousHistogram
    moments::Vector{T}
    values::Vector{T}
    errors::Vector{T}

    UniformHistogram{T}() where {T} = ( temp = @MArray zeros(T, 4); new( temp, zeros(T, 0), zeros(T, 0) ) )
    UniformHistogram(args...) = UniformHistogram{Float64}(args...)
end

KernelDistribution(::UniformHistogram) = UniformDistribution

"""
    uniform_function(::Number, μ, σ)
    uniform_function(::AbstractArray, μ, σ)

Calculate the normalized value of a uniform distribution centered at `val` that extends out by `err` above and below.
"""
uniform_function(x::Number, val, err) = ifelse( val - err < x && x < val + err, 1 / (2 * err), zero(val) )
uniform_function(x::AbstractArray, val, err) = broadcast( y -> uniform_function(y, val, err), x )
kernel(hist::UniformHistogram, x::AbstractArray, data) where {T} = uniform_function(x, data...) / length(hist)
