"""
    abstract type ContinuousHistogram end

`abstract type` representing the `ContinuousHistogram` concept. These are
histograms whose input data have _uncertainty_ associated with them, and therefore
we build them using a value-error-dependent kernel for each entry.

For each new `ContinuousHistogram`, one _must_ overload the functions given in each of the `src/Kernels` folder. Then one 
needs to add a new [`ContinuousDistribution`](@ref) and define its methods, like those shown
in the `src/Moments` files. Finally, the [`KernelDistribution`](@ref) must be mapped.
    
Then all `util`ity and `stats` functionality should _just work_.

!!! note
    In this context, _continuity_ refers to the domain of the histogram, and not necessarily its range.
"""
abstract type ContinuousHistogram end

"""
    KernelDistribution(::ContinuousHistogram) -> MethodError
    KernelDistribution(::GaussianHistogram)   -> GaussianDistribution
    KernelDistribution(::UniformHistogram)    -> UniformDistribution

Interface function that assigns a [`ContinuousHistogram`](@ref) `Type` to 
a specific kernel distribution.
"""
KernelDistribution(hist::ContinuousHistogram) = throw(MethodError(KernelDistribution, hist))

@doc raw"""
    kernel(::ContinuousHistogram, ::AbstractArray, data) -> MethodError
    kernel(::GaussianHistogram, ::AbstractArray, data)
    kernel(::UniformHistogram, ::AbstractArray, data)

The kernel function ``\mathcal{K}(x)`` used to compute a [`ContinuousHistogram`](@ref) such that

```math
\mathcal{H}(x) = \sum_{i = 1}^M \mathcal{K}_i(x),
```

where ``M`` is the total number of data points. For brevity, we have represented all `data`-dependence
through the subscript ``i``. 

!!! note
    By default, we `MethodError` out for any `<: ContinuousHistogram` until its directly implemented.
"""
kernel(hist::ContinuousHistogram, args...) = throw( MethodError(kernel, hist, args...) )

"""
    val_err(::ContinuousHistogram, idx)

Return a `(value, error)`-`Tuple`.
"""
val_err(hist::ContinuousHistogram, idx) = (hist.values[idx], hist.errors[idx])

"""
    construct!(output, ::ContinuousHistogram, x)

Similar to [`construct`](@ref) but here, the `output` `Array` is modified in-place.
"""
function construct!(output, hist::ContinuousHistogram, x)
    size(x) == size(output) ? nothing : ArgumentError("input and output vectors are of different sizes: $(size(x)) != $(size(output))")
    for idx ∈ 1:length(hist)
        @views output .+= kernel(hist, x, val_err(hist, idx))
    end
    return output
end

"""
    construct(::ContinuousHistogram, x)

Map the values of `x` through the [`ContinuousHistogram`] and return an
`Array` of the same `size` as `x`.
"""
function construct(hist::ContinuousHistogram, x)
    output = zeros(size(x))
    construct!(output, hist, x)
end

"""
    measurement(::ContinuousHistogram)

Interface to [`Measurements.jl`](https://juliaphysics.github.io/Measurements.jl/stable/). 
Return a `Measurement` with `val` as the [`ContinuousHistogram`](@ref) [`mean`](@ref) and the 
`err` as the [`ContinuousHistogram`](@ref) [`std`](@ref) (standard deviation).

!!! warning
    If the [`ContinuousHistogram`](@ref) in question is severely non-Gaussian, the first 
    two statistical [cumulants](https://en.wikipedia.org/wiki/Cumulant?oldformat=true) may
    be insufficient to appropriately describe the underlying distribution. In this sense, 
    a `measurement = value ± error` may not make sense to describe one's data.
"""
measurement(hist::ContinuousHistogram) = measurement(mean(hist), std(hist))