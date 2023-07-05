using Pkg
Pkg.activate(".")
Pkg.instantiate()

using Distributions
using DelimitedFiles
using ChromDistABC

const pp = ARGS[1] # data set


const pref=ARGS[2] # output prefix
const np = parse(Int, ARGS[3]) # n particles requested
const st = parse(Float64, ARGS[4])


const matAll = DelimitedFiles.readdlm(pp, ',')
const matY = matAll[matAll[:,3].!=1.0,:]
const matO = matAll[matAll[:,3].==1.0,:]
const countsY = map(x -> sum(matY[:,1].==x), 2:4)
const countsO = map(x -> sum(matO[:,1].==x), 2:4)
const totY = sum(countsY)
const totO=sum(countsO)
const distsY = matY[:,2]
const distsO = matO[:,2]

setupCounts = ABCRejection(simFunCounts, 
 		5, # number of parameters
 		12.0, #target ϵ
 		constants=DistConst((Y=countsY, O=countsO), st=st),
	nparticles=np,
 		Prior([Uniform(0.0, 1.0), # select model
           Uniform(0.1/20, 5.0/20), #α1young
		   Uniform(0.1/20, 10.0/20), #α2young
           Uniform(0.1/20, 5.0/20), #α1old
		   Uniform(0.1/20, 10.0/20)]); #α2old
 		maxiterations = 10^7, #Maximum number of iterations before the algorithm terminates 
 		)
 		
cOutSer = runabc(setupCounts, vcat(countsY, countsO), parallel=false)
writeoutput(cOutSer, file = pref*"outDist.abc")
