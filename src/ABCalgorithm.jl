#File which defined all the algorithms. Each algorithm takes in an ABCtype


"""
runabc(ABCsetup::ABCRejection, targetdata; progress = false, verbose = false, parallel = false)

Run ABC with ABCsetup defining the algorithm and inputs to algorithm, targetdata is the data we wish to fit the model to and will be used as an input for the simulation function defined in ABCsetup. If progress is set to `true` a progress meter will be shown. Inference will be run in parallel via multithreading if `parallel = true`. The environmental variable JULIA_NUM_THREADS needs to be set prior to launching a julia session.
"""
function runabc(ABCsetup::ABCRejection, targetdata; progress = false, verbose = false, parallel = false)

  #initalize array of particles
  particlesall = ParticleRejection[]
  
  if progress
    p = Progress(ABCsetup.nparticles, 1, "Running ABC rejection algorithm...", 30)
  end

  if parallel
    Printf.@printf("Preparing to run in parallel on %i processors\n", nthreads())

    particles = Array{ParticleRejection}(undef, ABCsetup.maxiterations)
    distvec = zeros(Float64, ABCsetup.maxiterations) #store distances in an array
    i = Atomic{Int64}(0)
    cntr = Atomic{Int64}(0)
    @threads for its = 1:ABCsetup.maxiterations

      if i[] > ABCsetup.nparticles
        break
      end

      #get new proposal parameters
      newparams = getproposal(ABCsetup.prior, ABCsetup.nparams)
      #simulate with new parameters
      dist, out = ABCsetup.simfunc(newparams, ABCsetup.constants, targetdata)
      #keep track of all particles incase we don't reach nparticles with dist < ϵ
      particlesall[its] = ParticleRejection(newparams, dist, out)

      #if simulated data is less than target tolerance accept particle
      if dist < ABCsetup.ϵ
        particles[its] = ParticleRejection(newparams, dist, out)
        distvec[its] = dist
        atomic_add!(i, 1)
      end
      atomic_add!(cntr,1)

    end
    # Remove particles that are still #undef and corresponding distances
    idx = [isassigned(particles,ii) for ii in eachindex(particles)]
    particles = particles[idx]
    distvec = distvec[idx]
    i = length(particles)    # Number of accepted particles
    its = cntr[]    # Total number of simulations

  else
    Printf.@printf("Preparing to run in serial on %i processor\n", 1)

    
    
    distvec = Float64[] #store distances in an array
    i = 1 #set particle indicator to 1
    its = 0 #keep track of number of iterations
    while (i < (ABCsetup.nparticles + 1)) & (its < ABCsetup.maxiterations)
      
      its += 1
      #get new proposal parameters
      newparams = getproposal(ABCsetup.prior, ABCsetup.nparams)
      
      #simulate with new parameters
      dist, out = ABCsetup.simfunc(newparams, ABCsetup.constants, targetdata)
      #keep track of all particles incase we don't reach nparticles with dist < ϵ
      #particlesall[its] = ParticleRejection(newparams, dist, out)
      
      #if simulated data is less than target tolerance accept particle
      
      if dist < ABCsetup.ϵ
        
        push!(particlesall, ParticleRejection(newparams, dist, out))
        push!(distvec, dist)
        i +=1
        if progress
          next!(p)
        end
      end
    end
    i -= 1    # Correct to total number of particels
  end
  
  

  if i < ABCsetup.nparticles
    @warn "Only accepted $(i) particles with ϵ < $(ABCsetup.ϵ). \n\tYou may want to increase ϵ or increase maxiterations. \n"
  end

  out = ABCrejectionresults(particlesall, its, ABCsetup, distvec)
  return out
end


function runabc(ABCsetup::ABCRejectionModel, targetdata; progress = false, verbose = false)

  ABCsetup.nmodels > 1 || error("Only 1 model specified, use ABCRejection method to estimate parameters for a single model")

  #initalize array of particles
  particles = Array{ParticleRejectionModel}(undef,
                ABCsetup.Models[1].nparticles)

  i = 1 #set particle indicator to 1
  its = 0 #keep track of number of iterations
  distvec = zeros(Float64, ABCsetup.Models[1].nparticles) #store distances in an array

  if progress
    p = Progress(ABCsetup.Models[1].nparticles, 1, "Running ABC rejection algorithm...", 30)
  end

  while (i < (ABCsetup.Models[1].nparticles + 1)) & (its < ABCsetup.Models[1].maxiterations)

    its += 1
    #sample uniformly from models
    model = rand(1:ABCsetup.nmodels)
    #get new proposal parameters
    newparams = getproposal(ABCsetup.Models[model].prior, ABCsetup.Models[model].nparams)
    #simulate with new parameters
    dist, out = ABCsetup.Models[model].simfunc(newparams, ABCsetup.Models[model].constants, targetdata)

    #if simulated data is less than target tolerance accept particle
    if dist < ABCsetup.Models[1].ϵ
      particles[i] = ParticleRejectionModel(newparams, model, dist, out)
      distvec[i] = dist
      i +=1
      if progress
        next!(p)
      end
    end
  end

  i > ABCsetup.Models[1].nparticles || error("Only accepted $(i) particles with ϵ < $(ABCsetup.Models[1].ϵ). \n\tDecrease ϵ or increase maxiterations ")

  out = ABCrejectionmodelresults(particles, its, ABCsetup, distvec)
  return out
end
