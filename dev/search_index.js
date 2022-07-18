var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = UncertainHistogramming","category":"page"},{"location":"#UncertainHistogramming","page":"Home","title":"UncertainHistogramming","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for UncertainHistogramming.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [UncertainHistogramming]","category":"page"},{"location":"#Base.show-Tuple{IO, GaussianHistogram}","page":"Home","title":"Base.show","text":"show([::IO,] ::GaussianHistogram)\nshow(::GaussianHistogram)\n\nprint the relevant information for a GaussianHistogram.\n\n\n\n\n\n","category":"method"},{"location":"#StatsBase.kurtosis","page":"Home","title":"StatsBase.kurtosis","text":"kurtosis(::GaussianHistogram [, excess = true ])\n\nPearson (excess) kurtosis of a GaussianHistogram.\n\nnote: Note\nFor comparison purposes, this kurtosis definition, with excess == true, applied to a  Gaussian distribution yields 0. If excess == false, then the kurtosis for a Gaussian is 3.\n\n\n\n\n\n","category":"function"},{"location":"#StatsBase.skewness-Tuple{GaussianHistogram}","page":"Home","title":"StatsBase.skewness","text":"skewness(::GaussianHistogram)\n\nFisher's skewness of a GaussianHistogram.\n\nnote: Note\nFor comparison purposes, this skewness definition applied to a Gaussian distribution yields identically zero.\n\n\n\n\n\n","category":"method"},{"location":"#UncertainHistogramming._online_mean-Tuple{Any, Any, Any}","page":"Home","title":"UncertainHistogramming._online_mean","text":"_online_mean(xnew, μold, nold)\n\nUpdate the old mean μold over nold elements with the new data point xnew. This function implements\n\nmu_n+1 = mu_n + fracx_n+1 - mu_nn+1\n\n\n\n\n\n","category":"method"},{"location":"#UncertainHistogramming._update_moment!-Tuple{GaussianHistogram, Any, Any, Any}","page":"Home","title":"UncertainHistogramming._update_moment!","text":"_update_moment!(::GaussianHistogram, moment_t, μ, σ)\n\nUpdate the GaussianHistogram's moment_t using an _online_mean with the inclusion of the value-error pair (μ, σ).\n\n\n\n\n\n","category":"method"},{"location":"#UncertainHistogramming._update_moments!-Tuple{GaussianHistogram, Any, Any}","page":"Home","title":"UncertainHistogramming._update_moments!","text":"_update_moments!(::GaussianHistogram, μ, σ)\n\nUpdate the the non-central moments of the GaussianHistogram online with the new value-error pair (μ, σ). This is done because each is non-central moment is the mean non-central moments of each GaussianDistribution.\n\n\n\n\n\n","category":"method"}]
}
