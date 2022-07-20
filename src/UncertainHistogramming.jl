module UncertainHistogramming

import Base: eltype, push!, length, getindex, setindex!, show
import StaticArrays: MArray, @MArray, @SArray
import Statistics: mean, var, std
import StatsBase: skewness, kurtosis
import Measurements: Measurement, measurement
using RecipesBase

export 
# Base overloads in ./util.jl
       push!, eltype, length, getindex, setindex!, show,
# Statistics and StatsBase overloads in ./stats.jl
       mean, var, std, skewness, kurtosis,
# Measurements.jl overloads
       measurement,
# ./Moments/Moments.jl overloads
       moment, FirstMoment, SecondMoment, ThirdMoment, FourthMoment,
# Generic UncertainHistogramming exports
       ContinuousHistogram, kernel, val_err, construct, construct!,
# Specific UncertainHistogramming types
       GaussianHistogram, UniformHistogram

include("Moments/Moments.jl")
include("ContinuousHistograms.jl")
include("util.jl")
include("stats.jl")
include("plotrecipes.jl")

# Specific kernel definitions (subtypes of ContinuousHistogram)
include("Kernels/GaussianKernel.jl")
include("Kernels/UniformKernel.jl")

end
