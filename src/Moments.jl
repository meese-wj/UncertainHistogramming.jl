
export moment, Moment, FirstMoment, SecondMoment, ThirdMoment, FourthMoment,
       ContinuousDistribution, GaussianDistribution, skewness, kurtosis

abstract type Moment end
struct FirstMoment <: Moment end
struct SecondMoment <: Moment end
struct ThirdMoment <: Moment end
struct FourthMoment <: Moment end

abstract type ContinuousDistribution end
struct GaussianDistribution <: ContinuousDistribution end

moment(d::Type{ContinuousDistribution}, m::Type{Moment}, params...) = MethodError(moment, d, m, params...)

moment(::Type{GaussianDistribution}, ::Type{FirstMoment},  μ, σ) = μ 
function moment(::Type{GaussianDistribution}, ::Type{SecondMoment}, μ, σ)
    μ, σ = promote(μ, σ)
    return μ^2 + σ^2
end
function moment(::Type{GaussianDistribution}, ::Type{ThirdMoment}, μ, σ)
    μ, σ = promote(μ, σ)
    return μ^3 + 3 * μ * σ^2
end
function moment(::Type{GaussianDistribution}, ::Type{FourthMoment}, μ, σ)
    μ, σ = promote(μ, σ)
    return μ^4 + 6 * μ^2 * σ^2 + 3 * σ^4
end

skewness(::Type{GaussianDistribution}, μ, σ) = zero(μ)

function kurtosis(::Type{GaussianDistribution}, μ, σ)
    μ, σ = promote(μ, σ)
    return convert(typeof(μ), 3)
end
