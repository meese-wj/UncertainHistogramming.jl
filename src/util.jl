import Base: eltype, push!, length, getindex, setindex!, show

eltype(::GaussianHistogram{T}) where {T} = T
length(hist::GaussianHistogram) = length(hist.values)

function push!( hist::GaussianHistogram, tup::Tuple{Number, Number} )
    _update_moments!(hist, tup...)
    push!(hist.values, tup[1])
    push!(hist.errors, tup[2])
    return hist
end
push!(hist::GaussianHistogram, meas::Measurement) = push!(hist, (meas.val, meas.err))

function push!(hist::GaussianHistogram, tup_vec::Vector{Tuple{Number, Number}})
    for tup ∈ tup_vec
        push!(hist, tup)
    end
    return hist
end
function push!(hist::GaussianHistogram, meas_vec::Vector{Measurement})
    for meas ∈ meas_vec
        push!(hist, meas)
    end
    return hist
end

for moment_t ∈ moment_list
    @eval Base.getindex(hist::GaussianHistogram, ::Type{$moment_t}) = hist.moments[ MomentIndex($moment_t) ]
    @eval Base.setindex!(hist::GaussianHistogram{T}, val::S, ::Type{$moment_t}) where {T, S} = hist.moments[ MomentIndex($moment_t) ] = convert(T, val)
    @eval moment(hist::GaussianHistogram, ::Type{$moment_t}) = hist[$moment_t]
end

"""
    show([::IO,] ::GaussianHistogram)
    show(::GaussianHistogram)

`print` the relevant information for a [`GaussianHistogram`](@ref).
"""
function show(io::IO, hist::GaussianHistogram)
    println(io, "GaussianHistogram{$(eltype(hist))}:")
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
show(hist::GaussianHistogram) = show(stdout, hist)