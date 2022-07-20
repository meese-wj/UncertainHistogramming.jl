
"""
    UniformDistribution <: ContinuousDistribution

Trait representing a uniform distribution.
"""
struct UniformDistribution <: ContinuousDistribution end

@doc raw"""
    moment(::Type{UniformDistribution}, FirstMoment, val, err)
    
Analytic expression for the no-central `FirstMoment` from a [`UniformDistribution`](@ref):

```math
M_1 = \int_{-\infty}^{\infty} {\rm d}y\, \mathcal{U}(y; x, \epsilon)\cdot y = x. 
```
"""
moment(::Type{UniformDistribution}, ::Type{FirstMoment},  val, err) = val 

@doc raw"""
    moment(::Type{UniformDistribution}, SecondMoment, val, err)
    
Analytic expression for the no-central `SecondMoment` from a [`UniformDistribution`](@ref):

```math
M_2 = \int_{-\infty}^{\infty} {\rm d}y\, \mathcal{U}(y; x, \epsilon)\cdot y^2 = x^2 + \frac{1}{3}\epsilon^2. 
```
"""
function moment(::Type{UniformDistribution}, ::Type{SecondMoment}, val, err)
    val, err = promote(val, err)
    return val^2 + (err^2) / 3
end

@doc raw"""
    moment(::Type{UniformDistribution}, ThirdMoment, val, err)
    
Analytic expression for the no-central `ThirdMoment` from a [`UniformDistribution`](@ref):

```math
M_3 = \int_{-\infty}^{\infty} {\rm d}y\, \mathcal{U}(y; x, \epsilon)\cdot y^3 = x\left(x^2 + \epsilon^2 \right). 
```
"""
function moment(::Type{UniformDistribution}, ::Type{ThirdMoment}, val, err)
    val, err = promote(val, err)
    return val * (val^2 + err^2)
end

@doc raw"""
    moment(::Type{UniformDistribution}, FourthMoment, val, err)
    
Analytic expression for the no-central `FourthMoment` from a [`UniformDistribution`](@ref):

```math
M_4 = \int_{-\infty}^{\infty} {\rm d}y\, \mathcal{U}(y; x, \epsilon)\cdot y^4 = x^4 + 2x^2\epsilon^2 + \frac{1}{5}\epsilon^4. 
```
"""
function moment(::Type{UniformDistribution}, ::Type{FourthMoment}, val, err)
    val, err = promote(val, err)
    return val^4 + 2 * val^2 * err^2 + (err^4) / 5
end
