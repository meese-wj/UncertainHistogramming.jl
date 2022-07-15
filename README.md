# UncertainHistogramming

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://meese-wj.github.io/UncertainHistogramming.jl/stable) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://meese-wj.github.io/UncertainHistogramming.jl/dev) [![Build Status](https://github.com/meese-wj/UncertainHistogramming.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/meese-wj/UncertainHistogramming.jl/actions/workflows/CI.yml?query=branch%3Amain)

A lightweight Julia package to plot a density function for a given set of values with known uncertainties.

This package provides similar functionality to what is expected from [_kernel density estimation_ (KDE)](https://en.wikipedia.org/wiki/Kernel_density_estimation?oldformat=true), but here the data errors/uncertainties which act as the kernel bandwidths are all, in principle, different.
