
function copyparticle(particle::ParticleRejection)
  return ParticleRejection(copy(particle.params), copy(particle.distance),
  [1])
end

function copyparticle(particle::ParticleRejectionModel)
  return ParticleRejectionModel(copy(particle.params), copy(particle.model),
  copy(particle.distance), [1])
end

"""
    ksdist(x::AbstractVector{T}, y::AbstractVector{S}) where {T <: Real, S <: Real}

Compute Kolmogorov-Smirnaov distance between `x` and `y`.
"""
function ksdist(x::AbstractVector{T}, y::AbstractVector{S}) where {T <: Real, S <: Real}

  #adapted from HypothesisTest.jl
  n_x, n_y = length(x), length(y)
  sort_idx = sortperm([x; y])
  pdf_diffs = [ones(n_x)/n_x; -ones(n_y)/n_y][sort_idx]
  cdf_diffs = cumsum(pdf_diffs)
  δp = maximum(cdf_diffs)
  δn = -minimum(cdf_diffs)
  δ = max(δp, δn)

  return δ
end

function show(io::IO, ABCresults::ABCrejectionresults)

  upperci = zeros(Float64, size(ABCresults.parameters, 2))
  lowerci = zeros(Float64, size(ABCresults.parameters, 2))
  parametermeans = zeros(Float64, size(ABCresults.parameters, 2))
  parametermedians = zeros(Float64, size(ABCresults.parameters, 2))

  for i in 1:size(ABCresults.parameters, 2)
    parametermeans[i] = mean(ABCresults.parameters[:, i])
    parametermedians[i] = median(ABCresults.parameters[:, i])
    (lowerci[i], upperci[i]) = quantile(ABCresults.parameters[:, i], [0.025,0.975])
  end

  @printf("Number of simulations: %.2e\n", ABCresults.numsims)
  @printf("Acceptance ratio: %.2e\n\n", ABCresults.accratio)

  print("Median (95% intervals):\n")
  for i in 1:length(parametermeans)
      @printf("Parameter %d: %.2f (%.2f,%.2f)\n", i, parametermedians[i], lowerci[i], upperci[i])
  end
end

function show(io::IO, ABCresults::ABCrejectionmodelresults)

  @printf("Number of simulations: %.2e\n", ABCresults.numsims)
  @printf("Acceptance ratio: %.2e\n\n", ABCresults.accratio)
  print("Model frequencies:\n")
  for j in 1:length(ABCresults.modelfreq)
    @printf("\tModel %d: %.2f\n", j, ABCresults.modelfreq[j])
  end

  print("\nParameters:\n\n")

  for j in 1:length(ABCresults.parameters)
    print("Model $j\n")

    upperci = zeros(Float64, size(ABCresults.parameters[j], 2))
    lowerci = zeros(Float64, size(ABCresults.parameters[j], 2))
    parametermeans = zeros(Float64, size(ABCresults.parameters[j], 2))
    parametermedians = zeros(Float64, size(ABCresults.parameters[j], 2))

    for i in 1:size(ABCresults.parameters[j], 2)
      parametermeans[i] = mean(ABCresults.parameters[j][:, i])
      parametermedians[i] = median(ABCresults.parameters[j][:, i])
      (lowerci[i], upperci[i]) = quantile(ABCresults.parameters[j][:, i], [0.025,0.975])
    end

    print("\tMedian (95% intervals):\n")
    for i in 1:length(parametermeans)
        @printf("\tParameter %d: %.2f (%.2f,%.2f)\n", i, parametermedians[i], lowerci[i], upperci[i])
    end
  end
end

function writeoutput(results::ABCrejectionresults; dir = "", file = "Rejection-output.txt")
  distance = map(x -> x.distance, results.particles)
  nparams =  size(results.parameters)[2]

  head = map(x -> "parameter$x\t", 1:nparams)
  append!(head, ["distance\n"])

  out = hcat(results.parameters, distance)

  f = open(joinpath(dir, file), "w")
  write(f, "## ABC Rejection algorithm\n")
  write(f, "## Number of simulations: $(results.numsims)\n")
  write(f, "## Acceptance ratio: $(round(results.accratio, sigdigits = 4))\n")
  write(f, head...)
  writedlm(f, out)
  close(f)

end

function writeoutput(results::ABCrejectionmodelresults; dir = "", file = "RejectionModel-output")

    for i in 1:length(results.modelfreq)
      prtcles = results.particles[map(x -> x.model, results.particles).==i]
      distance = map(x -> x.distance, prtcles)
      nparams =  size(results.parameters[i])[2]

      head = map(x -> "parameter$x\t", 1:nparams)
      append!(head, ["distance\n"])
      out = hcat(results.parameters[i], distance)

      f = open(joinpath(dir, "$(file)model$i.txt"), "w")
      write(f, "## ABC Rejection algorithm\n")
      write(f, "## Number of simulations: $(results.numsims)\n")
      write(f, "## Acceptance ratio: $(round(results.accratio, sigdigits = 4))\n")
      write(f, head...)
      writedlm(f, out)
      close(f)
  end

  #write model probabilities to file
  head = ["Model\t", "Probability\n"]
  f = open(joinpath(dir, "$(file)modelprobabilities.txt"), "w")
  write(f, head...)
  writedlm(f, hcat(collect(1:length(results.modelfreq)), results.modelfreq))
  close(f)

end


"""
	euc3(v::Vec3)

Euclidean distance from origin (0,0,0) for `Vec3` objects.
"""
euc3(v::Vec3) = sqrt(v.x^2 + v.y^2 + v.z^2)

"""
    normalize(v::Vec3)

Returns normalised version of `v` i.e. with length 1.
"""
function normalize(v::Vec3)
	l = sqrt(v.x^2 + v.y^2 + v.z^2)
	return Vec3(v.x/l, v.y/l,	v.z/l)
end



"""
    goodFocus(fo::Vec3,r)

A function to test whether a potential focus is within a spherical cell with ratius `r`.
"""
function goodFocus(fo::Vec3,r)
	euc3(fo) <= r
end


# Arithmetic for Ve3 objects
+(v::Vec3, u::Vec3)::Vec3 = Vec3(v.x + u.x, v.y + u.y, v.z + u.z)
*(v::Vec3, s::Real)::Vec3 = Vec3(v.x * s, v.y * s, v.z * s)
/(v::Vec3, s::Real)::Vec3 = Vec3(v.x / s, v.y / s, v.z / s)


"""
    mobileFociFast5(r=10, pars=[1.5/20, 20, 3.5/20, 20])

Generate four loci in two pairs inside a cell of radius `r`. The focus distance is sampled from scaled Beta2 distribution defined by `pars` (the slots of which are [μ₁, ν₁, μ₂, ν₂]). This implementation is using structs instead of vectors.
"""
function mobileFociFast5(r=10, pars=[1.5/20, 20, 3.5/20, 20])
	

	while true
		# 4 foci (Vec3)
		# 3 distances (Vec3)
		# 3 midpoints (Vec3)
		# 3 directions (scalars)

		f1 = rand3(-r, r)
		if goodFocus(f1, r)

	
		#distance to 2nd
		dist1 = rand(Beta2(pars[1], pars[2]))*2r
		
		# direction
		dirTemp = rand3(-1,1)
		dir1 = normalize(dirTemp)
		
		# 2nd locus
	
		f2 = f1 + (dir1 * dist1)
		if goodFocus(f2, r)
		mps1 = (f1 + f2) /2
		dist2 = rand(Beta2(pars[3], pars[4]))*2r
	    
		dirTmp = rand3(-1,1)
		dir2 = normalize(dirTmp)
		
		## 3rd locus
		mps2 = mps1 + (dir2 * dist2)
	#	
		dist3 = rand(Beta2(pars[1], pars[2]))*2r
	    dirTmp = rand3(-1, 1)
		dir3 = normalize(dirTmp)
		## normalise direction
		f3 = mps2 + (dir3 * dist3 / 2)
		if goodFocus(f3, r)
		## 4th locus
		f4 = mps2 + (dir3 * dist3 / -2)
		if goodFocus(f4, r)
			#return f1, f2, f3, f4
			return [f1.x f1.y f1.z; f2.x f2.y f2.z; f3.x f3.y f3.z; f4.x f4.y f4.z]
			end
		end #good focus 3
		end # if good focus 2
	end # if good focus 1
		

	end # while
end # function
	
"""
randFociFast5(r=10)

Generate four iid loci within a spherical cell of radius `r`.
"""
function randFociFast5(r=10)
	

	while true
		f1 = rand3(-r, r)
		if goodFocus(f1, r)
			f2 = rand3(-r, r)
			if goodFocus(f2, r)
				f3 = rand3(-r, r)
				if goodFocus(f3, r)
					f4 = rand3(-r, r)
					if goodFocus(f4, r)
						return [f1.x f1.y f1.z; f2.x f2.y f2.z; f3.x f3.y f3.z; f4.x f4.y f4.z]
					end
				end
			end
		end
	end # while
end # function
	
"""
    noDiscernible2(foci, minDist=2)

Given an m × 3 array of focus positions and a minimum distance, retuns the number of discernible signals. Before running hierarchical clustering, checks whether there are any pairwise distances < `minDist`.
"""
function noDiscernible2(foci, minDist=2)
	dd=pairwise(Euclidean(), foci')
	if any((dd[LowerTriangular(dd) .> 0.0] .< 0.5))
		ee=hclust(dd)
		return length(unique(cutree(ee; h=minDist)))
	else
		return size(foci)[1]
	end
	
end


"""
    getLociWithinSlice(foci, w=6, st=0)

Select those loci that are winthin a slice of width=`w` and with a top at `st`
"""
function getLociWithinSlice(foci, w=6, st=0)
	
	foci[ [(x[3] < st) & (x[3] > (st-w)) for x in eachrow(foci)], :]
end



function rand3(min, max)::Vec3
	Vec3(rand(Uniform(min, max)), rand(Uniform(min, max)), rand(Uniform(min, max)))
end

"""
    Beta2(μ, ν)

Re-parametrised version of `Distributions.Beta()`
"""
Beta2(μ, ν) = Beta(μ*ν, (1 − μ)*ν)


"""
    simFunCounts(params, constants, targetTab)

Simulation function for count data.
"""
function simFunCounts(params, constants, targetTab)
	
  # select model
  simY = Int[]
  simO = Int[]
  if params[1] < 1/4 # two params
    while length(simY) < constants[1]
      a = noDiscernible2(getLociWithinSlice(randFociFast5(10), 6, rand()*26-10), 0.5) 
      a >= 2 && push!(simY, a)
    end # while

    while length(simO) < constants[2]
      a = noDiscernible2(getLociWithinSlice(randFociFast5(10), 6, rand()*26-10), 0.5)
      a >= 2 && push!(simO, a)
    end

    return cityblock(targetTab, vcat(map(x->sum(simY.==x), 2:4), map(x->sum(simO.==x), 2:4))), 1

  elseif params[1] < 2/4	
      α1y=params[2]
      α2y=params[3]
      α1o=params[2]
      α2o=params[3]
  elseif params[1] < 3/4
      α1y=params[2]
      α2y=params[3]
      α1o=params[4]
      α2o=params[3]
  else
      α1y=params[2]
      α2y=params[3]
      α1o=params[4]
      α2o=params[5]
  end # if

  while length(simY) < constants[1]
    a = noDiscernible2(getLociWithinSlice(mobileFociFast5(10, [α1y, 20, α2y, 20]), 6, rand()*26-10), 0.5) 
    a >= 2 && push!(simY, a)
  end # while
  
  while length(simO) < constants[2]
    a = noDiscernible2(getLociWithinSlice(mobileFociFast5(10, [α1o, 20, α2o, 20]), 6, rand()*26-10), 0.5) 
    a >= 2 && push!(simO, a)
  end # while
  
  cityblock(targetTab, vcat(map(x->sum(simY.==x), 2:4), map(x->sum(simO.==x), 2:4))), 1
end


