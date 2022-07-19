
@recipe function f(x, hist::ContinuousHistogram; nkernels = 5)

    legend --> false
    linecolor --> :blue
    linewidth --> 2
    seriestype --> :path

    y = construct(hist, x)

    @series begin
        seriestype := :path
        linecolor := nothing
        fillcolor := :blue
        fillalpha := 0.075
        fillrange := 0.0
        x, y
    end

    if length(hist) <= nkernels
        for idx ∈ 1:min(length(hist), nkernels)
            @series begin
                seriestype := :path 
                linecolor := :orange 
                linestyle := :dash
                x, kernel(hist, x, val_err(hist, idx))
            end
        end
    end

    xlabel --> "\$x\$"
    ylabel --> "\$\\mathcal{H}(x)\$"
    
    x, y
end