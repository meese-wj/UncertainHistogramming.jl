# UncertainHistogramming 

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://meese-wj.github.io/UncertainHistogramming.jl/stable) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://meese-wj.github.io/UncertainHistogramming.jl/dev) [![Build Status](https://github.com/meese-wj/UncertainHistogramming.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/meese-wj/UncertainHistogramming.jl/actions/workflows/CI.yml?query=branch%3Amain)

A lightweight Julia package to plot a density function for a given set of values with known uncertainties.

An example application of the main `export`ed `abstract` `struct`, `ContinuousHistogram`, is to visualize a "histogram" of experimental values, when each _value_ has a measured experimental _uncertainty_. This is to be contrast with normal `Histogram`ming that assumes that each _value_ is exact, meaning its _uncertainty_ is `zero()`.

This package provides similar functionality to what is expected from [_kernel density estimation_ (KDE)](https://en.wikipedia.org/wiki/Kernel_density_estimation?oldformat=true), but here the data errors/uncertainties which act as the kernel bandwidths are all, in principle, different.

An example `GaussianHistogram <: ContinuousHistogram` can be `construct`ed from the following Julia code

```julia
using Plots # One must Pkg.add this separately
using UncertainHistogramming
Plots.gr()  # Use GR to reproduce the plot exactly 

# Create a list of (value, error)-tuples
values_errs = [(-3.5, 0.5), (-1.5, 0.75),
               (0, 0.25), (1.5, 0.75), 
               (3.5, 0.5)]

# Initialize the GaussianHistogram and 
# push! the data into it
hist = GaussianHistogram()
push!(hist, values_errs)

# Create an argument array to build the
# ContinuousHistogram from 
x = LinRange(-6, 6, 3000)

# Plot using a User Recipe 
# (thanks to RecipesBase.jl)
plot(x, hist)
```
![Example GaussianHistogram](docs/src/assets/example_gaussian_histogram.png)