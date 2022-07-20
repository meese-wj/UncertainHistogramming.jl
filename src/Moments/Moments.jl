
"""
    abstract type Moment end

Concept representing non-central statistical moments.
"""
abstract type Moment end

"""
    MomentIndex(::Type{<: Moment}) = -1

Mapping from `<:` [`Moment`](@ref) traits to indexable integers. 
"""
MomentIndex(::Type{<: Moment}) = -1

"""
    const moment_list::SVector

Convenience list of all moment type symbols used. Each of these 
`Symbol`s are `@eval`uated into an empty `struct` trait `<:` [`Moment`](@ref).

```@repl
UncertainHistogramming.moment_list
```
"""
const moment_list = @SArray [ :FirstMoment, :SecondMoment, :ThirdMoment, :FourthMoment ]
for (idx, moment_t) ∈ enumerate(moment_list)
    @eval struct $moment_t <: Moment end
    @eval MomentIndex(::Type{$moment_t}) = $idx
end

"""
    abstract type ContinuousDistribution end

Concept establishing what a _continuous distribution_ is within the code base.
Here, as with [`ContinuousHistogram`](@ref), the _continuity_ refers to the domain
of the distribution, and not necessarily its range.
"""
abstract type ContinuousDistribution end

moment(d::Type{ContinuousDistribution}, m::Type{Moment}, params...) = MethodError(moment, d, m, params...)

include("GaussianMoments.jl")
include("UniformMoments.jl")