function priorprob(parameters::Array{Float64, 1}, prior::Prior)
  pprob = 1
  for i in 1:length(parameters)
      pprob = pprob * pdf(prior.distribution[i], parameters[i])
  end
  return pprob
end

function kernelprob(p1, p2, kernel::Kernel)
    prob = 1
    for i in 1:length(p1.params)
      prob = prob * kernel.pdf_function(p1, p2, kernel.kernel_parameters, i)
    end
    return prob
end

function getmodelfreq(particles, ABCsetup)
  freq = zeros(Int64, ABCsetup.nmodels)
  models = map(x -> x.model, particles)
  for i in 1:ABCsetup.nmodels
    freq[i] = sum(models.==i)
  end
  return freq
end

function getmodelprob(currmodel, prevmodel, modelprob, ABCsetup)
  prob = ABCsetup.modelkern
  if currmodel == prevmodel
    return prob
  elseif sum(modelprob.>0.0) > 1
    return (1 - prob) / (sum(modelprob.>0.0) - 1)
  else
    return prob
  end
end


