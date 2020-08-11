
using Cxx
using Libdl

# to the path of ROPTLIB
# cd("/home/ubuntu/Documents/Sync/Codes/newROPTLIB/ROPTLIB")
cd("/Users/whuang/Documents/Syn/Codes/newROPTLIB/ROPTLIB")

# import the directory of the header files of Julia
# addHeaderDir("/home/ubuntu/Documents/julia-1.3.1/include/julia", kind=C_System) 
addHeaderDir("/Applications/Julia-1.3.app/Contents/Resources/julia/include/julia", kind=C_System) 

const path_to_lib = pwd()

# import the directories of ROPTLIB

addHeaderDir(path_to_lib * "/", kind=C_System) 
addHeaderDir(path_to_lib * "/BinaryFiles/", kind=C_System) 
addHeaderDir(path_to_lib * "/Manifolds/", kind=C_System) 
addHeaderDir(path_to_lib * "/Matlab/", kind=C_System) 
addHeaderDir(path_to_lib * "/Matlab/ForCpp/", kind=C_System) 
addHeaderDir(path_to_lib * "/Matlab/ForMatlab/", kind=C_System) 
addHeaderDir(path_to_lib * "/Others/", kind=C_System) 
addHeaderDir(path_to_lib * "/Problems/", kind=C_System) 
addHeaderDir(path_to_lib * "/Solvers/", kind=C_System) 
addHeaderDir(path_to_lib * "/test/", kind=C_System) 
addHeaderDir(path_to_lib * "/cwrapper/blas/", kind=C_System) 
addHeaderDir(path_to_lib * "/cwrapper/lapack/", kind=C_System) 

Libdl.dlopen(path_to_lib * "/DriverJuliaProb.so", Libdl.RTLD_GLOBAL)

cxx"""#define DRIVERJULIAPROB"""

cxxinclude("DriverJuliaProb.h")

# define the struct of parameters of solvers
# The meanings of the parameters can be found in Appendix B of the user manual

mutable struct SolverParams
	IsCheckParams::Int64
	IsCheckGradHess::Int64
	name::Ptr{UInt8}
	Stop_Criterion::Int64
	Tolerance::Float64
	Diffx::Float64
	NumExtraGF::Int64
	TimeBound::Int64
	Min_Iteration::Int64
	Max_Iteration::Int64
	OutputGap::Int64
	Verbose::Int64
	isconvex::Int64
	nu::Float64
	mu::Float64
	LengthSY::Int64
	lambdaLower::Float64
	lambdaUpper::Float64
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
mutable struct FunHandles
	fname::Ptr{UInt8}
	gfname::Ptr{UInt8}
	hfname::Ptr{UInt8}
	isstopped::Ptr{UInt8}
	LinesearchInput::Ptr{UInt8}
end

# define the struct of parameters of manifold
mutable struct ManiParams
	IsCheckParams::Int64
	numoftypes::Int64
	name::Ptr{UInt8}
	numofmani::Ptr{Int64}
	paramset::Ptr{Int64}
	m::Ptr{Int64}
	n::Ptr{Int64}
	p::Ptr{Int64}
end

# No default parameters of FunHandles and ManiParams are given.
# See JTestSimpleExample.jl for an example

# Call ROPTLIB API function
function DriverJuliaOPT(inHandles, inSparams, inMparams, inHasHHR, inX0)
    Sparamsvals = zeros(nfields(inSparams) - 1, 1);

    Sparamsvals[1] = inSparams.IsCheckParams
    Sparamsvals[2] = inSparams.IsCheckGradHess
    Sparamsvals[3] = inSparams.Stop_Criterion
    Sparamsvals[4] = inSparams.Tolerance
    Sparamsvals[5] = inSparams.Diffx
    Sparamsvals[6] = inSparams.NumExtraGF
    Sparamsvals[7] = inSparams.TimeBound
    Sparamsvals[8] = inSparams.Min_Iteration
    Sparamsvals[9] = inSparams.Max_Iteration
    Sparamsvals[10] = inSparams.OutputGap
    Sparamsvals[11] = inSparams.Verbose
    Sparamsvals[12] = inSparams.isconvex
    Sparamsvals[13] = inSparams.nu
    Sparamsvals[14] = inSparams.mu
    Sparamsvals[15] = inSparams.LengthSY
    Sparamsvals[16] = inSparams.lambdaLower
    Sparamsvals[17] = inSparams.lambdaUpper
    Sparamsvals[18] = inSparams.LineSearch_LS
    Sparamsvals[19] = inSparams.IsPureLSInput
    Sparamsvals[20] = inSparams.LS_alpha
    Sparamsvals[21] = inSparams.LS_beta
    Sparamsvals[22] = inSparams.Minstepsize
    Sparamsvals[23] = inSparams.Maxstepsize
    Sparamsvals[24] = inSparams.LS_ratio1
    Sparamsvals[25] = inSparams.LS_ratio2
    Sparamsvals[26] = inSparams.Initstepsize
    Sparamsvals[27] = inSparams.Accuracy
    Sparamsvals[28] = inSparams.Finalstepsize
    Sparamsvals[29] = inSparams.Num_pre_funs
    Sparamsvals[30] = inSparams.InitSteptype
    Sparamsvals[31] = inSparams.Acceptence_Rho
    Sparamsvals[32] = inSparams.Shrinked_tau
    Sparamsvals[33] = inSparams.Magnified_tau
    Sparamsvals[34] = inSparams.minimum_Delta
    Sparamsvals[35] = inSparams.maximum_Delta
    Sparamsvals[36] = inSparams.useRand
    Sparamsvals[37] = inSparams.Max_Inner_Iter
    Sparamsvals[38] = inSparams.Min_Inner_Iter
    Sparamsvals[39] = inSparams.theta
    Sparamsvals[40] = inSparams.kappa
    Sparamsvals[41] = inSparams.initial_Delta
    Sparamsvals[42] = inSparams.Eps
    Sparamsvals[43] = inSparams.Theta_eps
    Sparamsvals[44] = inSparams.Min_Eps
    Sparamsvals[45] = inSparams.Del
    Sparamsvals[46] = inSparams.Theta_del
	if(isreal(inX0))
        resultptr = @cxx DriverJuliaProb(inHandles.fname, inHandles.gfname, inHandles.hfname, inHandles.isstopped, inHandles.LinesearchInput,
                                         inSparams.name, pointer(Sparamsvals), nfields(inSparams) - 1,
                                         inMparams.name, inMparams.numoftypes, inMparams.numofmani, inMparams.paramset, inMparams.m, inMparams.n, inMparams.p, inMparams.IsCheckParams,
                                         inHasHHR, pointer(inX0), length(inX0))
	else
        resultptr = @cxx DriverJuliaProb(inHandles.fname, inHandles.gfname, inHandles.hfname, inHandles.isstopped, inHandles.LinesearchInput,
                                         inSparams.name, pointer(Sparamsvals), nfields(inSparams) - 1,
                                         inMparams.name, inMparams.numoftypes, inMparams.numofmani, inMparams.paramset, inMparams.m, inMparams.n, inMparams.p, inMparams.IsCheckParams,
                                         inHasHHR, convert(Ptr{Float64}, pointer(inX0)), length(inX0) * 2)
	end
	lentmp = unsafe_wrap(Array, resultptr, 1)
	resultlength = lentmp[1]
	
	resultArr = unsafe_wrap(Array, resultptr, convert(Int64, resultlength))
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
	lseries = convert(Int64, (resultlength - lx0 - 16) / 3)

	funs = view(resultArr, lx0 + 12 + 1 : lx0 + 12 + lseries)
	times = view(resultArr, lx0 + 12 + lseries + 1 : lx0 + 12 + 2 * lseries)
	grads = view(resultArr, lx0 + 12 + 2 * lseries + 1 : lx0 + 12 + 3 * lseries)

    Heigs = view(resultArr, lx0 + 12 + 3 * lseries + 1 : lx0 + 12 + 3 * lseries + 4)

	return (FinalIterate, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime,
			funs, grads, times, Heigs)
end

function real2complex(x::Array{Float64})::Array{Complex{Float64}}
	return unsafe_wrap(Array, convert(Ptr{Complex128}, pointer(x)), Int64(length(x) / 2), false)
end

function complex2real(x::Array{Complex{Float64}})::Array{Float64}
	return unsafe_wrap(Array, convert(Ptr{Float64}, pointer(x)), length(x) * 2, false)
end


