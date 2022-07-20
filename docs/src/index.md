```@meta
CurrentModule = UncertainHistogramming
```

# UncertainHistogramming [![GithubLink](assets/GitHub-Mark-Light-32px.png)](https://github.com/meese-wj/UncertainHistogramming.jl)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://meese-wj.github.io/UncertainHistogramming.jl/stable) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://meese-wj.github.io/UncertainHistogramming.jl/dev) [![Build Status](https://github.com/meese-wj/UncertainHistogramming.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/meese-wj/UncertainHistogramming.jl/actions/workflows/CI.yml?query=branch%3Amain)

Have you ever had a situation where you need to visualize a set of data as a histogram, except the data you have to visualize are each endowed with some amount of _uncertainty_? If so, this package is for you! [`UncertainHistogramming.jl`](https://github.com/meese-wj/UncertainHistogramming.jl) is a lightweight Julia package to plot a density function for a given set of values with known uncertainties.

## Background Information

An example application of the main `export`ed `abstract` `struct`, [`ContinuousHistogram`](@ref), is to visualize a "histogram" of experimental values, when each _value_ has a measured experimental _uncertainty_. This is to be contrast with normal `Histogram`ming that assumes that each _value_ is exact, meaning its _uncertainty_ is `zero()`.

For me, the need for this package first came about when I was running Monte Carlo simulations, where I needed to understand the underlying distribution of some observables. But, as anybody who has ever played around with Monte Carlo methods knows, each observable has a certain amount of statistical error. Thus, any regular histogram I would make when ignoring these statistical errors would not really expose the true distribution, as each data point could not entirely be claimed by a single histogram bin. So I invented the [`ContinuousHistogram`](@ref) as a somewhat tongue-in-cheek generalization of the regular histogram that takes data _uncertainty_ into account.

This package provides similar functionality to what is expected from [_kernel density estimation_ (KDE)](https://en.wikipedia.org/wiki/Kernel_density_estimation?oldformat=true), but here the data errors/uncertainties which act as the kernel bandwidths are all, in principle, different.

Currently, the only [`ContinuousHistogram`](@ref)/[`kernel`](@ref) pair are based on [`gaussian`](@ref)s, where each value-error pair are the mean and standard deviation of that [`gaussian`](@ref). In the future I plan to implement a couple of other subtypes of [`ContinuousHistogram`](@ref)s, but for right now, everything is _normal_. :-P

!!! note
    A [`ContinuousHistogram`](@ref) is _continuous_ in the sense of its _domain_. This is admittedly a bit confusion, but the discretization that occurs in a regular histogram comes from its _bins_, or its domain, _not_ its _range_. Of course, the range, or vertical values, are jumpy, but that is because of the discrete nature of the regular histogram. Most [kernel functions](https://en.wikipedia.org/wiki/Kernel_density_estimation?oldformat=true#Definition) that exist are at least piecewise continuous in their _range_, which is the same standard we take here.

## Example Usage

An example `GaussianHistogram <: ContinuousHistogram` can be `construct`ed from the following Julia code for a simple `Vector` of value-error `Tuple`s.

To start, first include the following packages:

```@example usage
using Plots # One must Pkg.add this separately
using UncertainHistogramming
Plots.gr()  # Use GR to reproduce the plot exactly
nothing # hide
```

Next, define a list of `(value, error)`-`Tuple`s:

```@example usage
values_errors = [(-3.5, 0.5), 
                 (-1.5, 0.75),
                 (0, 0.25), 
                 (1.5, 0.75), 
                 (3.5, 0.5)]
nothing # hide
```

From here, we're in the position to initialize a [`GaussianHistogram`](@ref) and [`push!`](@ref) the `values_errors` `Vector` into it.

```@example usage
hist = GaussianHistogram()
push!(hist, values_errors)
```

Note that the non-central statistical [`moment`](@ref)s are updated in an _online_ matter. This means that, aside from the overhead associated with [`push!`](@ref)ing two elements into the [`GaussianHistogram`](@ref)'s `values` and `errors` `Vectors`, there is an amortized cost associated with computing the statistics.

From here, we just need to define an input domain for the [`GaussianHistogram`](@ref) to be computed over as

```@example usage
x = LinRange(-6, 6, 3000)
```

and then, with the help of [`Plots.jl`](https://docs.juliaplots.org/stable/) and [`RecipesBase.jl`](https://juliaplots.org/RecipesBase.jl/stable/), we have

```@example usage
plot(x, hist)
```

In the plot, one can see that the [`GaussianHistogram`](@ref) is plotted as the solid blue curve, and the individual [`gaussian`](@ref) [`kernel`](@ref)s that make it up are plotted as the dashed orange curves.

!!! note
    One may also supply the keyword argument `nkernels` to `plot(x, hist)` to change the number of
    [`kernel`](@ref)s displayed. By default, `nkernels == 5`. 
    
    If the number of value-error pairs exceeds `nkernels`, that is `nkernels < length(hist)`, then no [`kernel`](@ref)s will be shown to save the end user from trying to understand and overly busy plot.



## Add [`UncertainHistogramming.jl`](https://github.com/meese-wj/UncertainHistogramming.jl) to your Julia environment

To add [`UncertainHistogramming.jl`](https://github.com/meese-wj/UncertainHistogramming.jl) simply press `]` in the Julia `REPL` to enter `pkg` mode and type

```
pkg> add UncertainHistogramming
```

and _presto!_ You now have full access to [`UncertainHistogramming.jl`](https://github.com/meese-wj/UncertainHistogramming.jl).
