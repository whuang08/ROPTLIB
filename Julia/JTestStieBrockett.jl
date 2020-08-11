
# set domain manifold to be the Stiefel manifold.
# Note that every parameter in ManiParams is an array. The idea to use an array
# is to make the framework compatible with produce of manifolds. See details in
# JTestProductExample.jl
mani = "Stiefel"; numoftypes = 1;
ms = [-1] # -1 means that m is not used in Stiefel manifold.
numofmani = [1]
# The size is R^{5 \times 3}
ns = [5];
ps = [3];
paramsets = [2];
Mparams = ManiParams(1, numoftypes, pointer(mani), pointer(numofmani), pointer(paramsets), pointer(ms), pointer(ns), pointer(ps))

# set function handles
fname = "func"
gfname = "gfunc"
hfname = "hfunc"
isstopped = "" # or "stopfunc" to use the function defined below
LinesearchInput = "" # or "LSfunc" to use the function defined below
Handles = FunHandles(pointer(fname), pointer(gfname), pointer(hfname), pointer(isstopped), pointer(LinesearchInput))

# set solvers by modifying the default one.
method = "RTRNewton"
Sparams.OutputGap = 1
Sparams.name = pointer(method)
Sparams.LineSearch_LS = 5
Sparams.Max_Iteration = 50
Sparams.IsPureLSInput = 0
Sparams.Stop_Criterion = 2
Sparams.IsCheckGradHess = 0
Sparams.IsCheckParams = 1

# use locking condition or not
HasHHR = 0

# Initial iterate and problem
import Random
Random.seed!(1)
n = ns[1]
p = ps[1]
B = randn(n, n)
B = B + B'
import SparseArrays
D = SparseArrays.spdiagm(0=>range(p, stop=1, length=p))
using LinearAlgebra
F = qr(randn(ns[1], ps[1]))
initialX = F.Q[:,1:p]

# Define function handles
# The function names are assigned to the "FunHandles" struct.
# See lines 17-21
function func(x, inTmp)
	x = reshape(x, n, p) # All the input argument is a vector. One has to reshape it to have a proper size
	outTmp = B * x * D
	fx = dot(x, outTmp)
	return (fx, outTmp) # The temparary data "outTmp" will replace the "inTmp"
end

function gfunc(x, inTmp) # The inTmp is the temparary data computed in "func".
	inTmp = reshape(inTmp, n, p) # All the input argument is a vector. One has to reshape it to have a proper size
	gf = 2.0::Float64 * inTmp
	return (gf, []) # If one does not want to change the temparary data, then let the outTmp be an empty array.
end

function hfunc(x, inTmp, eta)
	eta = reshape(eta, n, p) # All the input argument is a vector. One has to reshape it to have a proper size
	result = 2.0::Float64 * B * eta * D
	return (result, []) # If one does not want to change the temparary data, then let the outTmp be an empty array.
end

# Users can define their own stopping criterion by passing the name 
# of the function to the "isstopped" field in the object of structure FunHandles
function stopfunc(x, funs, ngfx, ngfx0)
# x: the current iterate
# gf: the gradient at x
# gx: the function value at x
# ngfx: the norm of gradient at x
# ngfx0: the norm of gradient at the initial iterate
	return (ngfx / ngfx0 < 1e-6)
end

# Users can define their own line search method by passing the name 
# of the function to the "LinesearchInput" field in the object of structure FunHandles
function LSfunc(x, eta, t0, s0)
# x: the current iterate
# eta: the search direction
# t0: the initial step size

# s0: the slope of the line search scalar function at zero
	return 1.0::Float64
end

# Call the solver and get results. See the user manual for details about the outputs.
#result = DriverJuliaOPT(Handles, Sparams, Mparams, HasHHR, initialX)
(FinalIterate, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times, Heigs) = DriverJuliaOPT(Handles, Sparams, Mparams, HasHHR, initialX);

