module ChromDistABC

export testFun


# Write your package code here.
"""
    testFun()

Prints "Test" to the terminal.
"""
testFun() = println("Test")




using Distributions
using ProgressMeter
using StatsBase
using RecipesBase
using Printf
using Distances
using DelimitedFiles
using Random
using Statistics
using Base.Threads

import Base.show
import LinearAlgebra: LowerTriangular
import Base.:+
import Base.:*
import Base.:/


export
  # types
  ABCtype,
  Prior,
  Particle,
  abctype,
  ParticleRejection,
  ABCRejection,
  ABCRejectionModel,
  Kernel,
  gaussiankernel,
  uniformkernel,

  #functions
  simFunCounts,
  ksdist,
  runabc,
  writeoutput

### source files
include("kernels.jl")
include("types.jl")
include("util.jl")
include("ABCalgorithm.jl")
include("sampling.jl")
include("calculateweights.jl")
include("plots.jl")


end
