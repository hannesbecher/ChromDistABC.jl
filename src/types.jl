abstract type ABCtype end
abstract type Particle end

"""
    Prior(distributions)

    Create Prior type for ABC algorithm specifying priors for each parameters. This is an array of Distribution types from Distribution.jl, each element corresponding to a parameter.
"""
mutable struct Prior
  distribution
  Prior(distributionarray) = new(tuple(distributionarray...))
end

mutable struct ParticleRejection <: Particle
  params::Array{Float64, 1}
  distance::Float64
  other::Any
end

mutable struct ParticleRejectionModel <: Particle
  params::Array{Float64, 1}
  model::Int64
  distance::Float64
  other::Any
end


"""
    ABCRejection(sim_func::Function, nparams::Int64, ϵ::Float64, prior::Prior; <keyword arguments>)

Create an ABCRejection type which will simulate data with sim_func. nparams is the number of parameters inputted into sim_func, ϵ is the target tolerance and prior sets the priors for the parameters. sim_func needs to take in 3 values, the parameters (in an array), constants (array) and target data in that order and needs to return 2 values, the first being the distance between the target data and simulated data and the second can be anything but is useful if for example you want to record some additional information about the simulations.
...
## Arguments
- `maxiterations = 10^5`: Maximum number of samples before the ABC algorithm terminates.
- `constants = []`: Any constants needed to simulate from sim_func
- `nparticles = 100`: Number of particles (ie samples) of ABC algorithm
...
"""
mutable struct ABCRejection <: ABCtype

  simfunc::Function
  nparams::Int64
  ϵ::Float64
  nparticles::Int64
  constants::Array{Any,1}
  maxiterations::Int64
  prior::Prior

  ABCRejection(sim_func::Function, nparams::Int64, ϵ::Float64, prior::Prior;
    maxiterations = 10000,
    constants = [],
    nparticles = 100,
    ) =
  new(sim_func, nparams, ϵ, nparticles, constants, maxiterations, prior)

end

"""
    ABCRejectionModel(sim_func::Array{Function, 1}, nparams::Array{Int64, 1}, ϵ::Float64, prior::Array{Prior, 1}; <keyword arguments>)

Create an ABCRejectionModel type which will create a type to run ABC with model selection. Each model is specified with a function, first input is an array of functions. nparams and priors are arrays for the number of parameters and priors for each model. each sim_func needs to take in 3 values, the parameters (in an array), constants (array) and target data in that order and needs to return 2 values, the first being the distance between the target data and simulated data and the second can be anything but is useful if for example you want to record some additional information about the simulations.
...
## Arguments
- `maxiterations = 10^5`: Maximum number of samples before the ABC algorithm terminates.
- `constants = [[]]`: Any constants needed to simulate from sim_func, needs to be an array of arrays, each one corresponding to a model function.
- `nparticles = 100`: Number of particles (ie samples) of ABC algorithm
...
"""
mutable struct ABCRejectionModel <: ABCtype

  Models::Array{ABCRejection, 1}
  nmodels::Int64

  ABCRejectionModel(sim_func::Array{Function, 1}, nparams::Array{Int64, 1}, ϵ::Float64, prior::Array{Prior, 1};
    constants = repeat([[]], outer = length(sim_func)),
    maxiterations = 10000,
    nparticles = 100,
    ) =
  new([ABCRejection(sim_func[i], nparams[i], ϵ, prior[i],  maxiterations = maxiterations, constants = constants[i], nparticles = nparticles) for i in 1:length(sim_func)], length(sim_func))

end

mutable struct ABCrejectionresults

  parameters::Array{Float64,2}
  accratio::Float64
  numsims::Int64
  dist::Array{Float64,1}
  particles::Array{ParticleRejection, 1}
  setup::ABCRejection

   function ABCrejectionresults(particles, its, ABCsetup, dist)
      parameters = hcat(map(x -> x.params, particles)...)'
      accratio = ABCsetup.nparticles/its
      new(parameters, accratio, its, dist, particles, ABCsetup)
   end
end

mutable struct ABCrejectionmodelresults

   parameters
   accratio::Float64
   numsims::Int64
   dist::Array{Float64,1}
   particles::Array{ParticleRejectionModel, 1}
   modelfreq::Array{Float64, 1}
   setup::ABCRejectionModel

   function ABCrejectionmodelresults(particles, its, ABCsetup, dist)
     parameters = []
     modelfreq = []
     for i in 1:ABCsetup.nmodels
        push!(modelfreq, sum(map(x -> x.model, particles) .== i))
        models = map(x -> x.model, particles)
        parameters = push!(parameters,
        hcat(map(x -> x.params, particles[map(x -> x.model, particles) .== i])...)')
     end
     accratio = ABCsetup.Models[1].nparticles/its
     modelfreq = modelfreq ./ sum(modelfreq)

     new(parameters, accratio, its, dist, particles, modelfreq, ABCsetup)
   end
end

"""
    Vec3(x, y, z)

Type for 3D coordinates.
"""
struct Vec3 
	x::Float64
	y::Float64
	z::Float64
end