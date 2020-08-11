
# set domain manifold to be the product of Stiefel manifolds: St(p, n)^2 \times St(q, m).
manis = "Stiefel, Stiefel"
numoftypes = 2
ms = [-1, -1] # -1 means that ms are not used in the Stiefel manifolds.
numofmani = [2, 1] # St(p, n) has power 2, therefore, the corresponding number is set to be 2.
# p = 3, n = 5, q = 2, m = 6
ns = [5, 6];
ps = [3, 2];
paramsets = [1, 1];
Mparams = ManiParams(1, numoftypes, pointer(manis), pointer(numofmani), pointer(paramsets), pointer(ms), pointer(ns), pointer(ps))

# set function handles
fname = "func_P"
gfname = "gfunc_P"
hfname = "hfunc_P"
isstopped = "stopfunc_P"
LinesearchInput = "LSfunc_P"
Handles = FunHandles(pointer(fname), pointer(gfname), pointer(hfname), pointer(isstopped), pointer(LinesearchInput))

# set solvers by modifying the default one.
method = "RTRNewton"
Sparams.name = pointer(method)
Sparams.OutputGap = 1
Sparams.LineSearch_LS = 5
Sparams.Max_Iteration = 50
Sparams.IsPureLSInput = 0
Sparams.IsCheckGradHess = 1

# use locking condition or not
HasHHR = 0

# Initial iterate and problem
import Random
Random.seed!(1)
n = ns[1]
p = ps[1]
m = ns[2]
q = ps[2]
B1 = randn(n, n)
B1 = B1 + B1'
import SparseArrays
D1 = SparseArrays.spdiagm(0=>range(p, stop=1, length=p))
B2 = randn(n, n)
B2 = B2 + B2'
D2 = SparseArrays.spdiagm(0=>range(p, stop=1, length=p))
B3 = randn(m, m)
B3 = B3 + B3'
D3 = SparseArrays.spdiagm(0=>range(q, stop=1, length=q))

using LinearAlgebra
F1 = qr(randn(n, p))
F2 = qr(randn(n, p))
F3 = qr(randn(m, q))
initialX1 = F1.Q[:,1:p]
initialX2 = F2.Q[:,1:p]
initialX3 = F3.Q[:,1:q]
initialX = [reshape(initialX1, n * p, 1); reshape(initialX2, n * p, 1); reshape(initialX3, m * q, 1)]

# Define function handles
# The function names are assigned to the "FunHandles" struct.
# See lines 17-21
function func_P(x, inTmp) # All the input argument is a vector.
	x1 = reshape(view(x, 1 : n * p), n, p)
	x2 = reshape(view(x, n * p + 1 : 2 * n * p), n, p)
	x3 = reshape(view(x, 2 * n * p + 1 : 2 * n * p + m * q), m, q)
	outTmp = [reshape(B1 * x1 * D1, n * p, 1); reshape(B2 * x2 * D2, n * p, 1); reshape(B3 * x3 * D3, m * q, 1)]
	fx = dot(x, outTmp)
	return (fx, outTmp) # The temparary data "outTmp" will replace the "inTmp"
end

function gfunc_P(x, inTmp)
	gf = 2.0::Float64 * inTmp
	return (gf, []) # If one does not want to change the temparary data, then let the outTmp be an empty array.
end

function hfunc_P(x, inTmp, eta)
	eta1 = reshape(view(eta, 1:n*p), n, p) # All the input argument is a vector. One has to reshape it to have a proper size
	eta2 = reshape(view(eta, n*p+1:2*n*p), n, p)
	eta3 = reshape(view(eta, 2*n*p+1:2*n*p+m*q), m, q)

	result = [reshape(2.0::Float64 * B1 * eta1 * D1, n * p, 1); reshape(2.0::Float64 * B2 * eta2 * D2, n * p, 1); reshape(2.0::Float64 * B3 * eta3 * D3, m * q, 1)]
	return (result, []) # If one does not want to change the temparary data, then let the outTmp be an empty array.
end

function stopfunc_P(x, funs, ngfx, ngfx0)
	return (ngfx / ngfx0 < 1e-6)
end

function LSfunc_P(x, eta, t0, s0)
	return 1.0::Float64
end

(FinalIterate, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times, Heigs) = DriverJuliaOPT(Handles, Sparams, Mparams, HasHHR, initialX)

