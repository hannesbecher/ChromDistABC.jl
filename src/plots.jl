
@recipe function f(results::ABCrejectionmodelresults)
    Probability = results.modelfreq
    Model = map(x -> "$x", 1:length(results.modelfreq))

    seriestype --> :bar
    title --> "Model Probabilities"
    yaxis --> "Probability"
    linecolor --> :white
    fillcolor --> :darkslategrey
    markerstrokecolor --> :white
    legend --> false
    grid --> false

    Model, Probability
end

@recipe function f(results::ABCrejectionmodelresults, model)
    nparams = size(results.parameters[model])[2]
    x = Array(results.parameters[model])
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
