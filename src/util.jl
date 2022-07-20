
"""
    eltype(::GaussianHistogram{T}) where {T} = T
    eltype(::UniformHistogram{T}) where {T} = T

`Base` overload for accessing the `eltype.`

!!! note
    One must overload this `util` function for all new [`ContinuousHistogram`] types.
"""
eltype(::GaussianHistogram{T}) where {T} = T
eltype(::UniformHistogram{T}) where {T} = T

"""
    length(::ContinuousHistogram)

`Base` overload to return the `length` of the `values` `Vector` in the argument [`ContinuousHistogram`](@ref).
"""
length(hist::ContinuousHistogram) = length(hist.values)

"""
    push!(::ContinuousHistogram, ::Tuple{Number, Number})
    push!(::ContinuousHistogram, ::Measurement)
    push!(::ContinuousHistogram, ::AbstractVector)
    push!(::ContinuousHistogram, ::Vector{Measurement})

`Base` overload to introduce a new value-error pair into the [`ContinuousHistogram`](@ref).
This function also [`_update_moments!`](@ref) in an amortized way to add little overhead.

!!! note
    The `AbstractVector` dispatch is really only meant for `Vector{Tuple{Number, Number}}`; 
    the latter of which contains no subtypes.
"""
function push!( hist::ContinuousHistogram, tup::Tuple{Number, Number} )
    _update_moments!(hist, tup...)
    push!(hist.values, tup[1])
    push!(hist.errors, tup[2])
    return hist
end
push!(hist::ContinuousHistogram, meas::Measurement) = push!(hist, (meas.val, meas.err))

function push!(hist::ContinuousHistogram, tup_vec::AbstractVector)
    for tup ∈ tup_vec
        push!(hist, tup)
    end
    return hist
end
function push!(hist::ContinuousHistogram, meas_vec::Vector{Measurement})
    for meas ∈ meas_vec
        push!(hist, meas)
    end
    return hist
end

@doc """
    getindex(::ContinuousHistogram, ::Type{<: Moment})

Convenience function to access a given [`Moment`](@ref) from its name.

```jldoctest
julia> hist = GaussianHistogram();

julia> push!(hist, (0, 1))
GaussianHistogram{Float64}:
  length  = 1
  moments = 0.0  1.0  0.0  3.0  

  Statistics
    mean        = 0.0
    variance    = 1.0
    skewness    = 0.0
    kurtosis    = 0.0

julia> hist[FirstMoment]
0.0

julia> hist[SecondMoment]
1.0

julia> hist[ThirdMoment]
0.0

julia> hist[FourthMoment]
3.0
```
"""
getindex

@doc """
    setindex!(::ContinuousHistogram, val, ::Type{<: Moment})

Convenience function for accessing the [`ContinuousHistogram`](@ref) [`Moment`](@ref)s from their name.
    
!!! warning
    This functionality should only be used internally as modifying the [`ContinuousHistogram`](@ref)
    `moments` directly would invalidate the [`push!`](@ref) pipeline and ultimately the statistics.
    
    But it's here if you need it for some `dev` reason.
"""
setindex!

@doc """
    moment(::ContinuousHistogram, ::Type{<: Moment})

Convenience wrapper to return a [`Moment`](@ref) from a [`ContinuousHistogram`](@ref) 
by that [`Moment`](@ref)'s name.

```jldoctest
julia> hist = GaussianHistogram();

julia> push!(hist, (0, 1))
GaussianHistogram{Float64}:
  length  = 1
  moments = 0.0  1.0  0.0  3.0  

  Statistics
    mean        = 0.0
    variance    = 1.0
    skewness    = 0.0
    kurtosis    = 0.0

julia> moment(hist, FirstMoment)
0.0

julia> moment(hist, SecondMoment)
1.0

julia> moment(hist, ThirdMoment)
0.0

julia> moment(hist, FourthMoment)
3.0
```
"""
moment

for moment_t ∈ moment_list
    @eval Base.getindex(hist::ContinuousHistogram, ::Type{$moment_t}) = hist.moments[ MomentIndex($moment_t) ]
    @eval Base.setindex!(hist::ContinuousHistogram, val, ::Type{$moment_t}) = hist.moments[ MomentIndex($moment_t) ] = convert(eltype(hist), val)
    @eval moment(hist::ContinuousHistogram, ::Type{$moment_t}) = hist[$moment_t]
end

"""
    show([::IO = stdout], ::ContinuousHistogram)
    show(::ContinuousHistogram)

`print` the relevant information for a [`ContinuousHistogram`](@ref).
"""
function show(io::IO, hist::ContinuousHistogram)
    println(io, "$(typeof(hist)):")
    println(io, "  length  = $(length(hist))")
    momstring = ""
    for mom ∈ hist.moments momstring *= "$mom  " end
    println(io, "  moments = $(momstring)")
    println(io, "")
    println(io, "  Statistics")
    println(io, "    mean        = $(mean(hist))")
    println(io, "    variance    = $(var(hist))")
    println(io, "    skewness    = $(skewness(hist))")
    println(io, "    kurtosis    = $(kurtosis(hist))")
end
show(hist::ContinuousHistogram) = show(stdout, hist)