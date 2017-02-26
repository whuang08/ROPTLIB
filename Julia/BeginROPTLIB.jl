
using Cxx

# to the path of ROPTLIB
cd("/home/whuang/Documents/ROPTLIB")

# import the directory of the header files of Julia
addHeaderDir("/home/whuang/Documents/julia/usr/include", kind=C_System) 
addHeaderDir("/home/whuang/Documents/julia/src", kind=C_System) 
addHeaderDir("/home/whuang/Documents/julia/src/support", kind=C_System) 

const path_to_lib = pwd()

# import the directories of ROPTLIB
addHeaderDir(path_to_lib * "/", kind=C_System) 
addHeaderDir(path_to_lib * "/BinaryFiles/", kind=C_System) 
addHeaderDir(path_to_lib * "/Manifolds/", kind=C_System) 
addHeaderDir(path_to_lib * "/Manifolds/CpxNStQOrth/", kind=C_System) 
addHeaderDir(path_to_lib * "/Manifolds/ElasticShape/", kind=C_System) 
addHeaderDir(path_to_lib * "/Manifolds/EucPositive/", kind=C_System) 
addHeaderDir(path_to_lib * "/Manifolds/Euclidean/", kind=C_System) 
addHeaderDir(path_to_lib * "/Manifolds/Grassmann/", kind=C_System) 
addHeaderDir(path_to_lib * "/Manifolds/L2Sphere/", kind=C_System) 
addHeaderDir(path_to_lib * "/Manifolds/LowRank/", kind=C_System) 
addHeaderDir(path_to_lib * "/Manifolds/Oblique/", kind=C_System) 
addHeaderDir(path_to_lib * "/Manifolds/OrthGroup/", kind=C_System) 
addHeaderDir(path_to_lib * "/Manifolds/PreShapeCurves/", kind=C_System) 
addHeaderDir(path_to_lib * "/Manifolds/SPDManifold/", kind=C_System) 
addHeaderDir(path_to_lib * "/Manifolds/SPDTensor/", kind=C_System) 
addHeaderDir(path_to_lib * "/Manifolds/Sphere/", kind=C_System) 
addHeaderDir(path_to_lib * "/Manifolds/Stiefel/", kind=C_System) 
addHeaderDir(path_to_lib * "/Matlab/", kind=C_System) 
addHeaderDir(path_to_lib * "/Matlab/ForCpp/", kind=C_System) 
addHeaderDir(path_to_lib * "/Matlab/ForCpp/Boundingbox/", kind=C_System) 
addHeaderDir(path_to_lib * "/Matlab/ForCpp/BrockettLRBFGSVTpaper/", kind=C_System) 
addHeaderDir(path_to_lib * "/Matlab/ForCpp/BrockettNonconvex/", kind=C_System) 
addHeaderDir(path_to_lib * "/Matlab/ForCpp/RepaRotCurves/", kind=C_System) 
addHeaderDir(path_to_lib * "/Matlab/ForCpp/SPDTensorDLandSC/", kind=C_System) 
addHeaderDir(path_to_lib * "/Matlab/ForCpp/SPDTensorDLandSC/DLandSC/", kind=C_System) 
addHeaderDir(path_to_lib * "/Matlab/ForCpp/SPDTensorDLandSC/EucPosSC/", kind=C_System) 
addHeaderDir(path_to_lib * "/Matlab/ForCpp/SPDTensorDLandSC/EucPosSC/algos/", kind=C_System) 
addHeaderDir(path_to_lib * "/Matlab/ForCpp/SPDTensorDLandSC/EucPosSC/tools/", kind=C_System) 
addHeaderDir(path_to_lib * "/Matlab/ForCpp/SPDTensorDLandSC/EucPosSC/tools/others/", kind=C_System) 
addHeaderDir(path_to_lib * "/Matlab/ForCpp/SPDTensorDLandSC/SPDtensorDL/", kind=C_System) 
addHeaderDir(path_to_lib * "/Matlab/ForCpp/SoftICA/", kind=C_System) 
addHeaderDir(path_to_lib * "/Matlab/ForCpp/SoftICA/New_folder/", kind=C_System) 
addHeaderDir(path_to_lib * "/Matlab/ForCpp/SoftICA/RBFGSNonconvexResults/", kind=C_System) 
addHeaderDir(path_to_lib * "/Matlab/ForCpp/SparsePCA/", kind=C_System) 
addHeaderDir(path_to_lib * "/Matlab/ForMatlab/", kind=C_System) 
addHeaderDir(path_to_lib * "/Matlab/ForMatlab/FromMelissa/", kind=C_System) 
addHeaderDir(path_to_lib * "/Others/", kind=C_System) 
addHeaderDir(path_to_lib * "/Problems/", kind=C_System) 
addHeaderDir(path_to_lib * "/Problems/ElasticCurvesRO/", kind=C_System) 
addHeaderDir(path_to_lib * "/Problems/EucFrechetMean/", kind=C_System) 
addHeaderDir(path_to_lib * "/Problems/EucPosSpCd/", kind=C_System) 
addHeaderDir(path_to_lib * "/Problems/EucQuadratic/", kind=C_System) 
addHeaderDir(path_to_lib * "/Problems/GrassRQ/", kind=C_System) 
addHeaderDir(path_to_lib * "/Problems/KarcherMean/", kind=C_System) 
addHeaderDir(path_to_lib * "/Problems/ObliqueTestSparsePCA/", kind=C_System) 
addHeaderDir(path_to_lib * "/Problems/OrthBoundingBox/", kind=C_System) 
addHeaderDir(path_to_lib * "/Problems/PreShapePathStraighten/", kind=C_System) 
addHeaderDir(path_to_lib * "/Problems/SPDMean/", kind=C_System) 
addHeaderDir(path_to_lib * "/Problems/SPDTensorDL/", kind=C_System) 
addHeaderDir(path_to_lib * "/Problems/ShapePathStraighten/", kind=C_System) 
addHeaderDir(path_to_lib * "/Problems/SphereConvexHull/", kind=C_System) 
addHeaderDir(path_to_lib * "/Problems/StieBrockett/", kind=C_System) 
addHeaderDir(path_to_lib * "/Problems/StieSoftICA/", kind=C_System) 
addHeaderDir(path_to_lib * "/Problems/StieSparseBrockett/", kind=C_System) 
addHeaderDir(path_to_lib * "/Problems/StieSumBrockett/", kind=C_System) 
addHeaderDir(path_to_lib * "/Problems/WeightedLowrank/", kind=C_System) 
addHeaderDir(path_to_lib * "/Solvers/", kind=C_System) 
addHeaderDir(path_to_lib * "/test/", kind=C_System) 
addHeaderDir(path_to_lib * "/cwrapper/blas/", kind=C_System) 
addHeaderDir(path_to_lib * "/cwrapper/lapack/", kind=C_System) 

Libdl.dlopen(path_to_lib * "/DriverJuliaProb.so", Libdl.RTLD_GLOBAL)

cxx"""#define DRIVERJULIAPROB"""

cxxinclude("DriverJuliaProb.h")

# define the struct of parameters of solvers
# The meanings of the parameters can be found in Appendix B of the user manual
type SolverParams
	IsCheckParams::Int64
	IsCheckGradHess::Int64
	name::Cstring
#Solvers
	Stop_Criterion::Int64
	Tolerance::Float64
	Diffx::Float64
	NumExtraGF::Int64
	TimeBound::Int64
	Min_Iteration::Int64
	Max_Iteration::Int64
	OutputGap::Int64
	DEBUG::Int64
#QuasiNewton
	isconvex::Int64
	nu::Float64
	mu::Float64
	LengthSY::Int64
	lambdaLower::Float64
	lambdaUpper::Float64
#SolversLS
	LineSearch_LS::Int64
	IsPureLSInput::Int64
	LS_alpha::Float64
	LS_beta::Float64
	Minstepsize::Float64
	Maxstepsize::Float64
	LS_ratio1::Float64
	LS_ratio2::Float64
	Initstepsize::Float64
	Accuracy::Float64
	Finalstepsize::Float64
	Num_pre_funs::Int64
	InitSteptype::Int64
#SolversTR
	Acceptence_Rho::Float64
	Shrinked_tau::Float64
	Magnified_tau::Float64
	minimum_Delta::Float64
	maximum_Delta::Float64
	useRand::Int64
	Max_Inner_Iter::Int64
	Min_Inner_Iter::Int64
	theta::Float64
	kappa::Float64
	initial_Delta::Float64
#SolversLSLPSub
	Eps::Float64
	Theta_eps::Float64
	Min_Eps::Float64
	Del::Float64
	Theta_del::Float64
end

# A default value for solver parameters is given:
method = "LRBFGS"
Sparams = SolverParams(1, 0, pointer(method),
		-1, -1, -1, -1, -1, -1, -1, -1, -1, # -1 means using the default parameters in the C++ code
		-1, -1, -1, -1, -1, -1,				# See the user manual for the setting details
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1)

# define the struct of parameters of function handles
# The strings are functions' names.
type FunHandles
	fname::Cstring
	gfname::Cstring
	hfname::Cstring
	isstopped::Cstring
	LinesearchInput::Cstring
end

# define the struct of parameters of manifold
type ManiParams
	IsCheckParams::Int64
	numoftypes::Int64
	name::Ptr{Cstring}
	numofmani::Ptr{Int64}
	paramset::Ptr{Int64}
	m::Ptr{Int64}
	n::Ptr{Int64}
	p::Ptr{Int64}
end

# No default parameters of FunHandles and ManiParams are given.
# See JTestSimpleExample.jl for an example

# Call ROPTLIB API function
function DriverJuliaOPT(inHandles, inSparams, inMparams, inHasHHR, inX0, inSoln = 0)
	if(isreal(inX0))
		ptrX0R::Ptr{Float64} = pointer(inX0)
        if(inSoln != 0)
            ptrSolnR::Ptr{Float64} = pointer(inSoln)
            resultptr = @cxx DriverJuliaProb(jpcpp"FunHandles"(inHandles), jpcpp"SolverParams"(inSparams), jpcpp"ManiParams"(inMparams), inHasHHR, ptrX0R, length(inX0), ptrSolnR)
        else
            resultptr = @cxx DriverJuliaProb(jpcpp"FunHandles"(inHandles), jpcpp"SolverParams"(inSparams), jpcpp"ManiParams"(inMparams), inHasHHR, ptrX0R, length(inX0))
        end
	else
		ptrX0C::Ptr{Float64} = convert(Ptr{Float64}, pointer(inX0))
        if(inSoln != 0)
            ptrSolnC::Ptr{Float64} = convert(Ptr{Float64}, pointer(inSoln))
            resultptr = @cxx DriverJuliaProb(jpcpp"FunHandles"(inHandles), jpcpp"SolverParams"(inSparams), jpcpp"ManiParams"(inMparams), inHasHHR, ptrX0C, length(inX0) * 2, ptrSolnC)
        else
            resultptr = @cxx DriverJuliaProb(jpcpp"FunHandles"(inHandles), jpcpp"SolverParams"(inSparams), jpcpp"ManiParams"(inMparams), inHasHHR, ptrX0C, length(inX0) * 2)
        end
	end
	lentmp = unsafe_wrap(Array, resultptr, 1, false)
	resultlength = lentmp[1]
	
	resultArr = unsafe_wrap(Array, resultptr, convert(Int64, resultlength), true)
	if(isreal(inX0))
		lx0 = length(inX0)
		Xopt = view(resultArr, 2:(lx0 + 1))
		FinalIterate = reshape(Xopt, size(inX0))
	else
		lx0 = length(inX0) * 2
		Xtmp::Array{Float64} = view(resultArr, 2:(lx0 + 1))
		Xopt = real2complex(Xtmp)
		FinalIterate = copy(reshape(Xopt, size(inX0)))
	end

	fv = resultArr[lx0 + 2]
	gfv = resultArr[lx0 + 3]
	gfgf0 = resultArr[lx0 + 4]
	iter = resultArr[lx0 + 5]
	nf = resultArr[lx0 + 6]
	ng = resultArr[lx0 + 7]
	nR = resultArr[lx0 + 8]
	nV = resultArr[lx0 + 9]
	nVp = resultArr[lx0 + 10]
	nH = resultArr[lx0 + 11]
	ComTime = resultArr[lx0 + 12]
	lseries = convert(Int64, (resultlength - lx0 - 12) / 4)

	funs = view(resultArr, lx0 + 12 + 1 : lx0 + 12 + lseries)
	grads = view(resultArr, lx0 + 12 + lseries + 1 : lx0 + 12 + 2 * lseries)
	times = view(resultArr, lx0 + 12 + 2 * lseries + 1 : lx0 + 12 + 3 * lseries)
	dists = view(resultArr, lx0 + 12 + 3 * lseries + 1 : lx0 + 12 + 4 * lseries)

	return (FinalIterate, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime,
			funs, grads, times, dists)
end

function real2complex(x::Array{Float64})::Array{Complex{Float64}}
	return unsafe_wrap(Array, convert(Ptr{Complex128}, pointer(x)), Int64(length(x) / 2), false)
end

function complex2real(x::Array{Complex{Float64}})::Array{Float64}
	return unsafe_wrap(Array, convert(Ptr{Float64}, pointer(x)), length(x) * 2, false)
end


