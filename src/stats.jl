mean( hist::GaussianHistogram ) = moment(hist, FirstMoment)
var( hist::GaussianHistogram ) = moment(hist, SecondMoment) - mean(hist)^2
std( hist::GaussianHistogram ) = sqrt( var(hist) )

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
