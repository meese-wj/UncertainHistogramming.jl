
"""
    mean(::ContinuousHistogram)

First moment of the [`ContinuousHistogram`](@ref).
"""
mean( hist::ContinuousHistogram ) = moment(hist, FirstMoment)
"""
    var(::ContinuousHistogram)

Second cumulant of the [`ContinuousHistogram`](@ref).
"""
var( hist::ContinuousHistogram ) = moment(hist, SecondMoment) - mean(hist)^2
"""
    std(::ContinuousHistogram)

The standard deviation of the [`ContinuousHistogram`](@ref).
"""
std( hist::ContinuousHistogram ) = sqrt( var(hist) )

"""
    skewness(::ContinuousHistogram)

[Fisher's skewness](https://en.wikipedia.org/wiki/Skewness?oldformat=true#Fisher's_moment_coefficient_of_skewness) of a [`ContinuousHistogram`](@ref).

!!! note
    For comparison purposes, this `skewness` definition applied to a Gaussian distribution yields identically zero.
"""
function skewness(hist::ContinuousHistogram)
    μ = mean(hist)
    val =  moment(hist, ThirdMoment) - 3 * μ * moment(hist, SecondMoment) + 2 * μ^3
    return val / std(hist)^3
end

"""
    kurtosis(::ContinuousHistogram [, excess = true ])

[Pearson (`excess`) kurtosis](https://en.wikipedia.org/wiki/Kurtosis?oldformat=true#Pearson_moments) of a [`ContinuousHistogram`](@ref).

!!! note
    For comparison purposes, this `kurtosis` definition, with `excess == true`, applied to a 
    Gaussian distribution yields `0`. If `excess == false`, then the `kurtosis` for a Gaussian is 3.
"""
function kurtosis(hist::ContinuousHistogram, excess = true)
    μ = mean(hist)
    val =  moment(hist, FourthMoment) - 4 * μ * moment(hist, ThirdMoment)
    val += 6 * μ^2 * moment(hist, SecondMoment) - 3 * μ^4
    val /= var(hist)^2
    return excess ? val - convert(typeof(val), 3) : val
end

"""
    _update_moments!(::ContinuousHistogram, val, err)

Update the the non-central moments of the [`ContinuousHistogram`](@ref) _online_ with
the new value-error pair `(val, err)`. This is done because each is non-central moment
is the mean non-central moments of each [`ContinuousDistribution`](@ref).
"""
function _update_moments!(hist::ContinuousHistogram, val, err)
    _update_moment!(hist, FirstMoment, val, err)
    _update_moment!(hist, SecondMoment, val, err)
    _update_moment!(hist, ThirdMoment, val, err)
    _update_moment!(hist, FourthMoment, val, err)
    return nothing
end

"""
    _update_moment!(::ContinuousHistogram, moment_t, val, err)

Update the [`ContinuousHistogram`](@ref)'s `moment_t` using an [`_online_mean`](@ref)
with the inclusion of the value-error pair `(val, err)`.
"""
@inline function _update_moment!(hist::ContinuousHistogram, moment_t, val, err) 
    hist[moment_t] = _online_mean( moment(KernelDistribution(hist), moment_t, val, err), 
                                   moment(hist, moment_t), length(hist) )
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
