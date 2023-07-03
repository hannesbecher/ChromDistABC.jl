module ChromDistABC

export testFun


# Write your package code here.
"""
    testFun()

Prints "Test" to the terminal.
"""
testFun() = println("Test")



using Clustering
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
import LinearAlgebra: LowerTriangular, diag
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

  #functions
  setupCounts,
  simFunCounts,
  simFunDists,
  ksdist,
  runabc,
  writeoutput

### source files
include("types.jl")
include("util.jl")
include("ABCalgorithm.jl")
include("sampling.jl")
include("plots.jl")


end
