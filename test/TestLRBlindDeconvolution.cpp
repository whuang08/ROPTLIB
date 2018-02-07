#include "test/TestLRBlindDeconvolution.h"

#ifdef ROPTLIB_WITH_FFTW

using namespace ROPTLIB;

//void testfft()
//{
//	int i, j, bw, bw2_1, size, size2_1, nrow, ncol;
//	int data_is_real;
//	int cutoff;
//	int rank, howmany_rank;
//	double *rresult, *iresult, *rdata, *idata;
//	double *workspace, *weights;
//
//	fftw_plan dctPlan;
//	fftw_plan fftPlan;
//	fftw_iodim dims[1], howmany_dims[1];
//
//	bw = 2;
//	weights = new double[4 * bw];// (double *)malloc(sizeof(double) * 4 * bw);
//	rdata = new double[5 * bw];// (double *)malloc(sizeof(double) * 5 * bw);
//	dctPlan = fftw_plan_r2r_1d(2 * bw, weights, rdata, FFTW_REDFT10, FFTW_ESTIMATE);
//	delete weights;
//	delete rdata;
//}

void testLRBlindDeconvolution(void)
{
	integer L = 20, K = 4, N = 5, r = 1;
	//integer m = 100, n = 15, r = 5;

	LowRank Domain(K, N, r);
	LowRankVariable InitialX(K, N, r);
	InitialX.RandInManifold();
	//Domain.CheckParams();

	//Domain.CheckIntrExtr(&InitialX);
	//Domain.CheckRetraction(&InitialX);
	//Domain.CheckcoTangentVector(&InitialX);
	//Domain.CheckDiffRetraction(&InitialX, false);
	//Domain.CheckIsometryofVectorTransport(&InitialX);

	//Domain.CheckLockingCondition(&InitialX);
	//Domain.CheckIsometryofInvVectorTransport(&InitialX);
	//Domain.CheckVecTranComposeInverseVecTran(&InitialX);
	//Domain.CheckTranHInvTran(&InitialX);
	//return;

	Euclidean EucDomain((K + N) * r);
	Domain.SetHasHHR(true);
	EucVariable EucX((K + N)*r);
	EucX.RandInManifold();

	//InitialX.Print("initialX:");

	// Generate the matrices in the Low rank approximation problem.
	double *y = new double[L * 2 + K * L * 2 + L * N * 2];
	double *B = y + L * 2;
	double *C = B + K * L * 2;
	for (integer i = 0; i < L * 2 + K * L * 2 + L * N * 2; i++)
		y[i] = genrandnormal();

	LRBlindDeconvolution Prob(y, B, 0, nullptr, nullptr, false, C, 0, nullptr, nullptr, 0, L, K, N, r);
	Prob.SetDomain(&Domain);

	EucBlindDeconvolution EucProb(y, B, 0, nullptr, nullptr, false, C, 0, nullptr, nullptr, 0, L, K, N, r);
	EucProb.SetDomain(&EucDomain);

	LRBFGS *RSDsolver = new LRBFGS(&EucProb, &EucX);
	RSDsolver->Debug = FINALRESULT;
	RSDsolver->OutputGap = 100;
	RSDsolver->Max_Iteration = 150;
	//RSDsolver->CheckParams();
	RSDsolver->Accuracy = 1e-6;
	RSDsolver->Finalstepsize = 1;
	RSDsolver->Tolerance = 1e-6;
	RSDsolver->Run();
	//EucProb.CheckGradHessian(&EucX);
	//EucProb.CheckGradHessian(RSDsolver->GetXopt());
	if (RSDsolver->Getnormgfgf0() < 1e-6)
		printf("SUCCESS!\n");
	else
		printf("FAIL!\n");
	delete RSDsolver;

	//RSDsolver = new LRBFGS(&Prob, &InitialX);
	////->LineSearch_LS = ARMIJO;
	////RSDsolver->LS_beta = 0.01;
	////RSDsolver->RCGmethod = DAI_YUAN;
	//RSDsolver->Debug = DETAILED;
	//RSDsolver->OutputGap = 100;
	//RSDsolver->Max_Iteration = 500;
	//RSDsolver->CheckParams();
	//RSDsolver->Accuracy = 1e-6;
	//RSDsolver->Finalstepsize = 1;
	//RSDsolver->Tolerance = 1e-10;
	//RSDsolver->Run();

	//// Compute the smallest eigenvalue of the Hessian at initial iterate.
	//SphereTx DomainPH0(Prob.GetDomain(), &InitialX);
	//SphereTxRQ ProbHess0(Prob.GetDomain(), &InitialX, &Prob, true);
	//ProbHess0.SetDomain(&DomainPH0);
	//Variable *TV00 = DomainPH0.RandominManifold();
	//RTRNewton *RTRNewtonsolver = new RTRNewton(&ProbHess0, TV00);
	//RTRNewtonsolver->Debug = NOOUTPUT;
	//RTRNewtonsolver->Run();
	//if (RTRNewtonsolver->Getnormgfgf0() > 1e-4)
	//	printf("Stop early when finding the smallest eigenvalue of the Hessian at initial iterate\n");
	//printf("Min eig at initial:%e\n", RTRNewtonsolver->Getfinalfun());
	//delete RTRNewtonsolver;
	//ProbHess0.SetMinorMax(false);
	//RTRNewtonsolver = new RTRNewton(&ProbHess0, TV00);
	//RTRNewtonsolver->Debug = NOOUTPUT;
	//RTRNewtonsolver->Run();
	//if (RTRNewtonsolver->Getnormgfgf0() > 1e-4)
	//	printf("Stop early when finding the largest eigenvalue of the Hessian at initial iterate\n");
	//delete TV00;
	//printf("Max eig at initial:%e\n", -RTRNewtonsolver->Getfinalfun());

	//// Compute the smallest eigenvalue of the Hessian at root.
	//Variable *root = InitialX.ConstructEmpty();
	//RSDsolver->GetXopt()->CopyTo(root);
	//SphereTx DomainPH(Prob.GetDomain(), root);
	//SphereTxRQ ProbHess(Prob.GetDomain(), root, &Prob, true);
	//ProbHess.SetDomain(&DomainPH);
	//Variable *TV0 = DomainPH.RandominManifold();
	//delete RTRNewtonsolver;
	//RTRNewtonsolver = new RTRNewton(&ProbHess, TV0);
	//RTRNewtonsolver->Debug = NOOUTPUT;
	//RTRNewtonsolver->Run();
	//if (RTRNewtonsolver->Getnormgfgf0() > 1e-4)
	//	printf("Stop early when finding the smallest eigenvalue of the Hessian at optimum\n");
	//printf("Min eig at Final:%e\n", RTRNewtonsolver->Getfinalfun());
	//delete RTRNewtonsolver;
	//ProbHess.SetMinorMax(false);
	//RTRNewtonsolver = new RTRNewton(&ProbHess, TV0);
	//RTRNewtonsolver->Debug = NOOUTPUT;
	//RTRNewtonsolver->Run();
	//if (RTRNewtonsolver->Getnormgfgf0() > 1e-4)
	//	printf("Stop early when finding the largest eigenvalue of the Hessian at optimum\n");
	//printf("Max eig at Final:%e\n", -RTRNewtonsolver->Getfinalfun());
	//delete RTRNewtonsolver;
	//delete root;
	//delete TV0;

	//delete RSDsolver;

	//RTRNewton *solver2 = new RTRNewton(&Prob, &InitialX);
	//solver2->Debug = ITERRESULT;
	//solver2->OutputGap = 10;
	//solver2->Tolerance = 1e-8;
	//solver2->Run();
	//delete solver2;

	delete[] y;
};

void testLRBlindDeconvolutionSparse(void)
{
	integer L = 20, K = 4, N = 5, r = 1;
	//integer m = 100, n = 15, r = 5;
	LowRank Domain(K, N, r);
	Domain.SetHasHHR(true);
	LowRankVariable InitialX(K, N, r);
	InitialX.RandInManifold();
	//Domain.CheckParams();
	//Domain.CheckIntrExtr(&InitialX);
	//Domain.CheckRetraction(&InitialX);
	//Domain.CheckcoTangentVector(&InitialX);
	//Domain.CheckDiffRetraction(&InitialX, false);
	//Domain.CheckIsometryofVectorTransport(&InitialX);

	//Domain.CheckLockingCondition(&InitialX);
	//Domain.CheckIsometryofInvVectorTransport(&InitialX);
	//Domain.CheckVecTranComposeInverseVecTran(&InitialX);
	//Domain.CheckTranHInvTran(&InitialX);
	//return;

	//InitialX.Print("initialX:");

	// Generate the matrices in the Low rank approximation problem.
	double *y = new double[L * 2 + K * L * 2 + L * N * 2];
	double *B = y + L * 2;
	double *C = B + K * L * 2;
	for (integer i = 0; i < L * 2 + K * L * 2 + L * N * 2; i++)
		y[i] = genrandnormal();
	integer nzmaxB = 2 * L * K;
	integer nzmaxC = 2 * L * N;
	int *irB = new int[4 * L * K + 4 * L * N];
	int *jcB = irB + 2 * L * K;
	int *irC = jcB + 2 * L * K;
	int *jcC = irC + 2 * L * N;
	for (integer i = 0; i < K; i++)
	{
		for (integer j = 0; j < 2 * L; j++)
		{
			irB[j + i * 2 * L] = j;
			jcB[j + i * 2 * L] = i;
		}
	}
	for (integer i = 0; i < N; i++)
	{
		for (integer j = 0; j < 2 * L; j++)
		{
			irC[j + i * 2 * L] = j;
			jcC[j + i * 2 * L] = i;
		}
	}
	Stiefel mani1(K, r);
	Euclidean mani2(r, r);
	Stiefel mani3(N, r);

	//blas_sparse_matrix sB;
	//sB = BLAS_duscr_begin(2 * L, K);
	////BLAS_duscr_insert_entries(sB, nzmaxB, B, irB, jcB);
	////BLAS_duscr_end(sB);
	//BLAS_usds(sB);


	LRBlindDeconvolution Prob(y, B, nzmaxB, irB, jcB, true, C, nzmaxC, irC, jcC, true, L, K, N, r);
	Prob.SetDomain(&Domain);

	//SphereTx DomainPH(&Domain, &InitialX);
	//SphereTxRQ ProbHess(&Domain, &InitialX, &Prob);
	//ProbHess.SetDomain(&DomainPH);
	//Vector *gf0 = DomainPH.RandominManifold();
	//ProbHess.CheckGradHessian(gf0);
	//delete gf0;


	//Prob.f(&InitialX);//---
	//InitialX.Print("X:");//---
	//ForDebug::Print("y:", y, 2 * L);//---
	//ForDebug::Print("B:", B, 2 * L, K);//--
	//ForDebug::Print("C:", C, 2 * L, N);//--

	//Prob.CheckGradHessian(&InitialX);

	LRBFGS *RSDsolver = new LRBFGS(&Prob, &InitialX);
	//->LineSearch_LS = ARMIJO;
	//RSDsolver->LS_beta = 0.01;
	//RSDsolver->RCGmethod = DAI_YUAN;
	RSDsolver->Debug = FINALRESULT;
	RSDsolver->OutputGap = 1;
	RSDsolver->Max_Iteration = 50;
	//RSDsolver->CheckParams();
	RSDsolver->Accuracy = 1e-6;
	RSDsolver->Finalstepsize = 1;
	RSDsolver->Tolerance = 1e-6;
	RSDsolver->Run();
	if (RSDsolver->Getnormgfgf0() < 1e-6)
		printf("SUCCESS!\n");
	else
		printf("FAIL!\n");
	//Prob.CheckGradHessian(&InitialX);//--
	//Prob.CheckGradHessian(RSDsolver->GetXopt());//--

	//// Compute the smallest eigenvalue of the Hessian at root.
	//Variable *root = InitialX.ConstructEmpty();
	//RSDsolver->GetXopt()->CopyTo(root);
	//SphereTx DomainPH(&Domain, root);
	//SphereTxRQ ProbHess(&Domain, root, &Prob, true);
	//ProbHess.SetDomain(&DomainPH);
	//Variable *TV0 = DomainPH.RandominManifold();
	//RTRNewton RTRNewtonsolver(&ProbHess, TV0);
	//RTRNewtonsolver.Debug = FINALRESULT;
	//RTRNewtonsolver.Run();
	//delete root;
	//delete TV0;

	//// Compute the smallest eigenvalue of the Hessian at initial iterate.
	//SphereTx DomainPH0(Prob.GetDomain(), &InitialX);
	//SphereTxRQ ProbHess0(Prob.GetDomain(), &InitialX, &Prob, true);
	//ProbHess0.SetDomain(&DomainPH0);
	//Variable *TV00 = DomainPH0.RandominManifold();
	//RTRNewton *RTRNewtonsolver = new RTRNewton(&ProbHess0, TV00);
	//RTRNewtonsolver->Debug = NOOUTPUT;
	//RTRNewtonsolver->Run();
	//if (RTRNewtonsolver->Getnormgfgf0() > 1e-4)
	//	printf("Stop early when finding the smallest eigenvalue of the Hessian at initial iterate\n");

	//delete RTRNewtonsolver;
	//ProbHess0.SetMinorMax(false);
	//RTRNewtonsolver = new RTRNewton(&ProbHess0, TV00);
	//RTRNewtonsolver->Debug = NOOUTPUT;
	//RTRNewtonsolver->Run();
	//if (RTRNewtonsolver->Getnormgfgf0() > 1e-4)
	//	printf("Stop early when finding the largest eigenvalue of the Hessian at initial iterate\n");
	//delete TV00;

	//// Compute the smallest eigenvalue of the Hessian at root.
	//Variable *root = InitialX.ConstructEmpty();
	//RSDsolver->GetXopt()->CopyTo(root);
	//SphereTx DomainPH(Prob.GetDomain(), root);
	//SphereTxRQ ProbHess(Prob.GetDomain(), root, &Prob, true);
	//ProbHess.SetDomain(&DomainPH);
	//Variable *TV0 = DomainPH.RandominManifold();
	//delete RTRNewtonsolver;
	//RTRNewtonsolver = new RTRNewton(&ProbHess, TV0);
	//RTRNewtonsolver->Debug = NOOUTPUT;
	//RTRNewtonsolver->Run();
	//if (RTRNewtonsolver->Getnormgfgf0() > 1e-4)
	//	printf("Stop early when finding the smallest eigenvalue of the Hessian at optimum\n");

	//delete RTRNewtonsolver;
	//ProbHess.SetMinorMax(false);
	//RTRNewtonsolver = new RTRNewton(&ProbHess, TV0);

	//RTRNewtonsolver->Debug = NOOUTPUT;
	//RTRNewtonsolver->Run();
	//if (RTRNewtonsolver->Getnormgfgf0() > 1e-4)
	//	printf("Stop early when finding the largest eigenvalue of the Hessian at optimum\n");
	//delete RTRNewtonsolver;
	//delete root;
	//delete TV0;

	delete RSDsolver;

	//RTRNewton *solver2 = new RTRNewton(&Prob, &InitialX);
	//solver2->Debug = ITERRESULT;
	//solver2->OutputGap = 10;
	//solver2->Tolerance = 1e-8;
	//solver2->Run();
	//delete solver2;

	delete[] y;
	delete[] irB;

};

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 7)
	{
		mexErrMsgTxt("The number of arguments should be at least seven.\n");
	}
	double *y, *B, *C, *X, *Xopt;
	y = mxGetPr(prhs[0]);
	B = mxGetPr(prhs[1]);
	C = mxGetPr(prhs[2]);
	X = mxGetPr(prhs[3]);
	/* dimensions of input matrices */
	integer L, K, N, HasHHR, r;
	L = mxGetM(prhs[0]) / 2;
	K = mxGetN(prhs[1]);
	N = mxGetN(prhs[2]);
	if (mxGetM(prhs[1]) != 2 * L || mxGetM(prhs[2]) != 2 * L)
	{
		mexErrMsgTxt("The size of B or the size of C is not correct.\n");
	}
	r = static_cast<integer> (mxGetScalar(prhs[4]));
	HasHHR = static_cast<integer> (mxGetScalar(prhs[5]));
	if (mxGetM(prhs[3]) != (K + N + r) * r)
	{
		mexErrMsgTxt("The size of initial x is not correct.\n");
	}

	bool isBsparse = mxIsSparse(prhs[1]);
	bool isCsparse = mxIsSparse(prhs[2]);
	integer nzmaxB = 0, nzmaxC = 0;
	size_t *irB = nullptr, *jcB = nullptr, *irC = nullptr, *jcC = nullptr;
	int *inirB = nullptr, *injcB = nullptr, *inirC = nullptr, *injcC = nullptr;
	if (isBsparse)
	{
		nzmaxB = mxGetNzmax(prhs[1]);
		irB = mxGetIr(prhs[1]);
		jcB = mxGetJc(prhs[1]);

		inirB = new int[2 * nzmaxB];
		injcB = inirB + nzmaxB;
		for (integer i = 0; i < K; i++)
		{
			for (unsigned long long j = jcB[i]; j < jcB[i + 1]; j++)
			{
				/*row: ir[j], column: i, entry: A[j]*/
				inirB[j] = irB[j];
				injcB[j] = i;
			}
		}
	}
	if (isCsparse)
	{
		nzmaxC = mxGetNzmax(prhs[2]);
		irC = mxGetIr(prhs[2]);
		jcC = mxGetJc(prhs[2]);

		inirC = new int[2 * nzmaxC];
		injcC = inirC + nzmaxC;
		for (integer i = 0; i < N; i++)
		{
			for (unsigned long long j = jcC[i]; j < jcC[i + 1]; j++)
			{
				/*row: ir[j], column: i, entry: A[j]*/
				inirC[j] = irC[j];
				injcC[j] = i;
			}
		}
	}

	genrandseed(0);
	CheckMemoryDeleted = new std::map<integer *, integer>;

	// Obtain an initial iterate by taking the Q factor of qr decomposition
	LowRankVariable LRX(K, N, r);
	double *LRXptr = LRX.ObtainWriteEntireData();
	for (integer i = 0; i < (K + N + r) * r; i++)
		LRXptr[i] = X[i];

	LowRankVariable *LRsoln = nullptr;
	if (nrhs >= 8)
	{
		double *soln = mxGetPr(prhs[7]); /*soln: n by p*/
		LRsoln = new LowRankVariable(K, N, r);
		double *LRsolnptr = LRsoln->ObtainWriteEntireData();
		for (integer i = 0; i < (K + N + r) * r; i++)
		{
			LRsolnptr[i] = soln[i];
		}
	}

	// Define the manifold
	LowRank Domain(K, N, r);

	// Define the matrix completion problem
	LRBlindDeconvolution Prob(y, B, nzmaxB, inirB, injcB, isBsparse, C, nzmaxC, inirC, injcC, isCsparse, L, K, N, r);
	Prob.SetDomain(&Domain);

	Domain.SetHasHHR(HasHHR != 0);
	//Domain.CheckParams();

	// Call the function defined in DriverMexProb.h
	ParseSolverParamsAndOptimizing(prhs[6], &Prob, &LRX, LRsoln, plhs);

	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	delete LRsoln;
	delete[] inirB;
	delete[] inirC;
	return;
}

#endif
#endif
