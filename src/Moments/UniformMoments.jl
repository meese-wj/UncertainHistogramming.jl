"""
    GaussianDistribution <: ContinuousDistribution

Trait representing a Gaussian distribution.
"""
struct GaussianDistribution <: ContinuousDistribution end

@doc raw"""
    moment(::Type{GaussianDistribution}, FirstMoment, μ, σ)
    
Analytic expression for the no-central `FirstMoment` from a [`GaussianDistribution`](@ref):

```math
M_1 = \int_{-\infty}^{\infty} {\rm d}y\, G(y;\mu,\sigma)\cdot y = \mu. 
```
"""
moment(::Type{GaussianDistribution}, ::Type{FirstMoment},  μ, σ) = μ 

@doc raw"""
    moment(::Type{GaussianDistribution}, SecondMoment, μ, σ)
    
Analytic expression for the no-central `SecondMoment` from a [`GaussianDistribution`](@ref):

```math
M_2 = \int_{-\infty}^{\infty} {\rm d}y\, G(y;\mu,\sigma)\cdot y^2 = \mu^2 + \sigma^2. 
```
"""
function moment(::Type{GaussianDistribution}, ::Type{SecondMoment}, μ, σ)
    μ, σ = promote(μ, σ)
    return μ^2 + σ^2
end

@doc raw"""
    moment(::Type{GaussianDistribution}, ThirdMoment, μ, σ)
    
Analytic expression for the no-central `ThirdMoment` from a [`GaussianDistribution`](@ref):

```math
M_3 = \int_{-\infty}^{\infty} {\rm d}y\, G(y;\mu,\sigma)\cdot y^3 = \mu^3 + 3\mu\sigma^2. 
```
"""
function moment(::Type{GaussianDistribution}, ::Type{ThirdMoment}, μ, σ)
    μ, σ = promote(μ, σ)
    return μ^3 + 3 * μ * σ^2
end

@doc raw"""
    moment(::Type{GaussianDistribution}, FourthMoment, μ, σ)
    
Analytic expression for the no-central `FourthMoment` from a [`GaussianDistribution`](@ref):

```math
M_4 = \int_{-\infty}^{\infty} {\rm d}y\, G(y;\mu,\sigma)\cdot y^4 = \mu^4 + 6\mu^2\sigma^2 + 3\sigma^4. 
```
"""
function moment(::Type{GaussianDistribution}, ::Type{FourthMoment}, μ, σ)
    μ, σ = promote(μ, σ)
    return μ^4 + 6 * μ^2 * σ^2 + 3 * σ^4
end
