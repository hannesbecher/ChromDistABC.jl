




@recipe function f(results::ABCrejectionresults)
    nparams = size(results.parameters)[2]
    x = Array(results.parameters)
    seriestype --> :histogram
    nbins --> 20
    layout --> nparams
    title --> hcat(map(x -> "Parameter $x", 1:nparams)...)
    linecolor --> :white
    fillcolor --> :darkred
    markerstrokecolor --> :white
    legend --> false
    x
end
