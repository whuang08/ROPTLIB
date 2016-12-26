
# JTestPhaseRetrieval.jl
# Apply Riemannian method for Phase Retrieval problem
#
# --- by Wen Huang

function JCompareInitialPhaseRetrieval()
	global n1 = 64
	global n2 = 64
	global n = n1 * n2
	global l = 5
	global m = l * n
	global p = 1
	global delta = 0.95
	global kappa = 1 / n
	global tol = 1e-9

	# problem related parameters
	global seed = floor(rand() * 1000000)
	#seed = 673358
	println("seed:", seed)
	srand(Int64(seed))
	global masks = complex(randn(n1, n2, l), randn(n1, n2, l))
#	global masks = octanaryrandom(n1, n2, l)
	global phase = complex(randn(n1, n2), randn(n1, n2))
	phase = phase / norm(phase)
	global vphase = reshape(phase, n, 1)
	global b = zeros(m, 1)
	for i in 1 : l
		zi = reshape(fft(phase .* masks[:, :, i]) / sqrt(n), n, 1)
		b[(i - 1) * n + 1 : i * n] = conj(zi) .* zi;
	end

	# initial iterate for Riemannian approach.
	(initialX, t) = InitialIterate(p, 0)
	Xrot = initialX * (initialX'*vphase) / norm(initialX'*vphase)
	println("0: dist(sol, xinit) / |sol|:", vecnorm(vphase - Xrot) / vecnorm(phase))
	# initial iterate for Riemannian approach.
	(initialX, t) = InitialIterate(p, Int64(floor(50 / p)))
	println("50: dist(sol, xinit) / |sol|:", vecnorm(vphase - initialX * (initialX'*vphase) / norm(initialX'*vphase)) / vecnorm(phase))
	
	println("\n")
	# initial iterate for Riemannian approach.
	p = 2; WFiter = 7; Rieiter = 7
	(initialX, inittime) = InitialIterate(p, WFiter)

	# initial step size at the first iterate
	(Xopt, nfftrank2, time) = JTestPhaseRetrievalfixedrank(initialX, Rieiter, WFiter * 2 * l * p)
	println("Rie:dist(sol, xinit) / |sol|:", vecnorm(vphase - Xopt * (Xopt'*vphase) / norm(Xopt'*vphase)) / vecnorm(phase))
	println("total time:", inittime + time)

	println("\n")
	# initial iterate for Riemannian approach.
	p = 1; WFiter = 16; Rieiter = 16
	(initialX, inittime) = InitialIterate(p, WFiter)

	# initial step size at the first iterate
	(Xopt, nfftrank2, time) = JTestPhaseRetrievalfixedrank(initialX, Rieiter, WFiter * 2 * l * p)
	println("Rie:dist(sol, xinit) / |sol|:", vecnorm(vphase - Xopt * (Xopt'*vphase) / norm(Xopt'*vphase)) / vecnorm(phase))
	println("total time:", inittime + time)
end

function octanaryrandom(n1, n2, l)
	masks = complex(zeros(n1, n2, l), zeros(n1, n2, l))
	for i in 1 : n1
		for j in 1 : n2
			for k in 1 : l
				tmp = rand()
				if(tmp < 0.25)
					b1 = 1
				elseif(tmp < 0.5)
					b1 = -1
				elseif(tmp < 0.75)
					b1 = -im
				else
					b1 = im
				end
				tmp = rand()
				if(tmp < 0.8)
					b2 = sqrt(2.0)/2
				else
					b2 = sqrt(3.0)
				end
				masks[i, j, k] = b1 * b2
			end
		end
	end
	return masks
end

function InitialIterate(k, niter)
	srand(Int64(seed))
	tic()
	lambda = sqrt(sum(b) / real(sum(conj(masks) .* masks)) * n)
#	initialX = qr(complex(randn(n, k), randn(n, k)))[1]
#	for i in 1 : niter
#		initialX = qr(Ax(initialX, k))[1]
#	end
	initialX = complex(randn(n, k), randn(n, k))
	initialX = initialX / norm(initialX)
	for i in 1 : niter
		initialX = Ax(initialX, k)
	end
	if(niter >= 1)
		initialX = initialX / norm(initialX)
	end
	time = toc()
	return (lambda * initialX, time)
end

function WFlow(initialX, tol, initnfft)
	gradf = egf(initialX)
	ngf = norm(gradf)
	ngf0 = ngf
	t = 0
	deno = norm(initialX)^2
	X = initialX
	while(ngf > tol && t < 350)
		# this stepsize works well for n = 2^2 to 256^2
		stepsize = min(1.0 - exp(-(1.0 + t) / 330.0), 0.036 * log(n)-0.015) * (exp(2.35*n^(1/8)-1.7)) / deno
		gradf = egf(X)
		if(mod(t, 10) == 0)
			ngf = norm(gradf)
			println("iter", t, ",ngf:", ngf, ", ngf/ngf0:", ngf / ngf0, ", stepsize:", stepsize)
		end
		X -= stepsize * gradf
		t+=1
	end
	ngf = norm(gradf)
	println("final iter:", t, ",ngf:", ngf, ", ngf/ngf0:", ngf / ngf0, ", nfft:", initnfft + (t+1)*2*l)
	return (X, initnfft + (t+1)*2*l)
end

function Ax(x, k)
	sqrtn = sqrt(n)
	ZY = complex(zeros(m, k), zeros(m, k))
	for i in 1 : l
		for j in 1 : k
			temp = reshape(fft(reshape(x[:, j], n1, n2) / sqrtn .* masks[:, :, i]), n, 1)
			ZY[(i - 1) * n + 1 : i * n, j] = temp
		end
	end
	DZY = complex(zeros(m, k), zeros(m, k))
	for i in 1 : k
		DZY[:, i] = b .* ZY[:, i]
	end

	Y = complex(zeros(n, k), zeros(n, k))
	for i in 1 : l
		temp = complex(zeros(n, k), zeros(n, k))
		for j in 1 : k
			temp[:, j] = reshape(ifft(reshape(DZY[(i - 1) * n + 1 : i * n, j] * sqrtn, n1, n2)), n, 1)
			temp[:, j] = conj(reshape(masks[:, :, i], n, 1)) .* temp[:, j]
		end
		Y += temp
	end
	return Y
end

#--------------------main test function for Riemannian method-----------------------------------------
function JTestPhaseRetrievalfixedrank(initialX, maxiter, initnfft)
	# set domain manifold to be the complex noncompact Stiefel manifold quotient the unitary group
	mani1 = "CpxNStQOrth"; ManArr = [pointer(mani1)]
	UseDefaultArr = [-1] # -1 means that the default value in C++ is used.
	numofmani = [1]
	# The size is n by p
	ns = [n]
	ps = [p]
	Mparams = ManiParams(0, length(ManArr), pointer(ManArr), pointer(numofmani), pointer(UseDefaultArr), pointer(UseDefaultArr), pointer(ns), pointer(ps))

	# set function handles
	fname = "func"
	gfname = "gfunc"
	hfname = "hfunc"
	isstopped = "stopfunc"
	LinesearchInput = "" # or "LSfunc" to use the function defined below
	Handles = FunHandles(pointer(fname), pointer(gfname), pointer(hfname), pointer(isstopped), pointer(LinesearchInput))

	# set solvers by modifying the default one.
	method = "LRBFGS"
	Sparams.name = pointer(method)
	Sparams.OutputGap = 10
	Sparams.LengthSY = 2
	Sparams.Min_Iteration = maxiter;
	Sparams.Max_Iteration = maxiter;
	Sparams.Tolerance = tol
	Sparams.Stop_Criterion = 2
	Sparams.DEBUG = 1
	Sparams.Accuracy = 1e-6
	Sparams.Finalstepsize = 1
	Sparams.IsCheckGradHess = 0
	Sparams.IsCheckParams = 0
	Sparams.isconvex = 1
	#Sparams.Initstepsize = stepsize0

	# use locking condition or not
	HasHHR = 0::Int64
	nfft = initnfft; nn = 0;
	# Call the solver and get results. See the user manual for details about the outputs.
	(FinalIterate, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times) = DriverJuliaOPT(Handles, Sparams, Mparams, HasHHR, initialX)

	nfft = (nf + ng) * l * p; nn = 4 * ng * p^2 + 4 * (length(funs) - 1) * p^2 + nV * 5 * p^2;
	(U, S) = svd(FinalIterate)

	return (U[:, 1] * S[1], nfft, ComTime)
end


# Define function handles
# The function names are assigned to the "FunHandles" struct.
# See lines 17-21
function func(xreal, inTmpreal)
	x = real2complex(xreal)
	x = reshape(x, n, p) # All the input argument is a vector. One has to reshape it to have a proper size
	sqrtn = sqrt(n)
	outTmp = complex(zeros(m, p + 1), zeros(m, p + 1)) # the first k columns are for ZY and the last column is for D.
	for i in 1 : l
		sZqi = zeros(n, 1)
		for j in 1 : p
			temp = reshape(fft(reshape(x[:, j], n1, n2) / sqrtn .* masks[:, :, i]), n, 1)
			outTmp[(i - 1) * n + 1 : i * n, j] = temp
			sZqi += conj(temp) .* temp
		end
		outTmp[(i - 1) * n + 1 : i * n, p + 1] = sZqi - b[(i - 1) * n + 1 : i * n]
	end
	output = norm(outTmp[:, p + 1])^2 / norm(b)^2
#	if(p != 1)
#		output += kappa * vecnorm(x)^2
#	end
	return (output, complex2real(outTmp)) # The temparary data "outTmp" will replace the "inTmp"
end

function gfunc(xreal, inTmpreal) # The inTmp is the temparary data computed in "func".
	x = real2complex(xreal)
	inTmp = real2complex(inTmpreal)
	x = reshape(x, n, p) # All the input argument is a vector. One has to reshape it to have a proper size
	inTmp = reshape(inTmp, m, p + 1) # All the input argument is a vector. One has to reshape it to have a proper size
	sqrtn = sqrt(n)
	DZY = complex(zeros(m, p), zeros(m, p))
	for i = 1 : p
		DZY[:, i] = inTmp[:, end] .* inTmp[:, i]
	end

	Y = complex(zeros(n, p), zeros(n, p))
	for i in 1 : l
		temp = complex(zeros(n, p), zeros(n, p))
		for j in 1 : p
			temp[:, j] = reshape(ifft(reshape(DZY[(i - 1) * n + 1 : i * n, j] * sqrtn, n1, n2)), n, 1)
			temp[:, j] = conj(reshape(masks[:, :, i], n, 1)) .* temp[:, j]
		end
		Y += temp
	end
	Y = 4.0 * Y / norm(b)^2
#	if(p != 1)
#		Y += kappa * 2 * x
#	end
	gf = Y / (x' * x)

	return (complex2real(gf), []) # If one does not want to change the temparary data, then let the outTmp be an empty array.
end

# action of the Hessian is not used. So use identity to avoid error.
function hfunc(xreal, inTmpreal, etareal)
	return (etareal, []) # If one does not want to change the temparary data, then let the outTmp be an empty array.
end

# Users can define their own stopping criterion by passing the name 
# of the function to the "isstopped" field in the object of structure FunHandles
function stopfunc(xreal, gfreal, fx, ngfx, ngfx0)
# x: the current iterate
# gf: the gradient at x
# gx: the function value at x
# ngfx: the norm of gradient at x
# ngfx0: the norm of gradient at the initial iterate
	return false
end


function egf(x) # the Euclidean gradient for Wirtinger flow method
	sqrtn = sqrt(n)
	ZY = complex(zeros(m, p), zeros(m, p))
	D = zeros(m, 1)
	for i in 1 : l
		sZqi = zeros(n, 1)
		for j in 1 : p
			temp = reshape(fft(reshape(x[:, j], n1, n2) / sqrtn .* masks[:, :, i]), n, 1)
			ZY[(i - 1) * n + 1 : i * n, j] = temp
			conjtemp = conj(temp);
			sZqi += conj(temp) .* temp
		end
		D[(i - 1) * n + 1 : i * n] = sZqi - b[(i - 1) * n + 1 : i * n]
	end

	DZY = complex(zeros(m, p), zeros(m, p))
	for i = 1 : p
		DZY[:, i] = D .* ZY[:, i]
	end
	Y = complex(zeros(n, p), zeros(n, p))
	for i in 1 : l
		temp = complex(zeros(n, p), zeros(n, p))
		for j in 1 : p
			temp[:, j] = reshape(ifft(reshape(DZY[(i - 1) * n + 1 : i * n, j] * sqrtn, n1, n2)), n, 1)
			temp[:, j] = conj(reshape(masks[:, :, i], n, 1)) .* temp[:, j]
		end
		Y += temp
	end
	Y = 4.0 * Y / norm(b)^2
	return Y
end

