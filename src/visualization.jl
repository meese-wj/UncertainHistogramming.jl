
# Create an empty struct for dispatch
# which only has the field args and 
# lowercase plotting function and seriestype
@userplot UncertainHistogram 

@recipe function f(uh::UncertainHistogram)
    if length(uh.args) != 2 || !(typeof(h.args[1]) <: AbstractVector) || !(typeof(h.args[2]) <: ContinuousHistogram)
        throw(ArgumentError("Uncertain Histograms must be given a Vector and ContinuousHistogram. Got: $(typeof(uh.args))"))
    end

    # Extract out the args
    # construct the ContinuousHistogram
    x, hist = uh.args
    y = construct(hist, x)

    legend --> false
    linecolor --> :blue
    seriestype --> :path

    # Main plot
    @series begin
        x, y
    end
end