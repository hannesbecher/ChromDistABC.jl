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

    
    i = Atomic{Int64}(0) # counting up to target number of particles
    cntr = Atomic{Int64}(0) # counting overall loops

    @threads for its = 1:ABCsetup.maxiterations

      if i[] > ABCsetup.nparticles
        break
      end

      #get new proposal parameters
      newparams = getproposal(ABCsetup.prior, ABCsetup.nparams)
      #simulate with new parameters
      dist, out = ABCsetup.simfunc(newparams, ABCsetup.constants, targetdata)
      #keep track of all particles incase we don't reach nparticles with dist < ϵ
      

      #if simulated data is less than target tolerance accept particle
      if dist < ABCsetup.ϵ
        push!(particlesall, ParticleRejection(newparams, dist, out))
        atomic_add!(i, 1)
      end
      atomic_add!(cntr,1)

    end
    
    
    
    
    i = length(particlesall)    # Number of accepted particles
    its = cntr[]    # Total number of simulations

  else
    Printf.@printf("Preparing to run in serial on %i processor\n", 1)

    
    
    
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

  distvec = map(x -> x.distance, particlesall)
  out = ABCrejectionresults(particlesall, its, ABCsetup, distvec)
  return out
end

