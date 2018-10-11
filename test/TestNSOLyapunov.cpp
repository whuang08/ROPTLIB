#include "test/TestNSOLyapunov.h"

using namespace ROPTLIB;

void testNSOLyapunov(void)
{
	integer n = 50, p = 2;

	double *A = new double[n * n * 3];
	double *M = A + n * n;
	double *C = M + n * n;

	double *Temp = new double[n * n];
	for (integer i = 0; i < n * n; i++)
		Temp[i] = genrandnormal();
	// A = temp * temp^T, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
	dgemm_(GLOBAL::N, GLOBAL::T, &n, &n, &n, &GLOBAL::DONE, Temp, &n, Temp, &n, &GLOBAL::DZERO, A, &n);

	for (integer i = 0; i < n * n; i++)
		M[i] = 0;
	for (integer i = 0; i < n; i++)
		M[i + i * n] = 1;
	//for (integer i = 0; i < n * n; i++)
	//	Temp[i] = genrandnormal();
	//// M = temp * temp^T, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
	//dgemm_(GLOBAL::N, GLOBAL::T, &n, &n, &n, &GLOBAL::DONE, Temp, &n, Temp, &n, &GLOBAL::DZERO, M, &n);

	for (integer i = 0; i < n * n; i++)
		Temp[i] = genrandnormal();
	// C = temp * temp^T, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
	dgemm_(GLOBAL::N, GLOBAL::T, &n, &n, &GLOBAL::IONE, &GLOBAL::DONE, Temp, &n, Temp, &n, &GLOBAL::DZERO, C, &n);

	delete[] Temp;

	NSOVariable SCOX(n, p);

	SCOX.RandInManifold();

	// Define the manifold
	NStQOrth Domain(n, p);
	//Domain.SetHasHHR(true);

	Domain.ChooseNSOParamsSet1();

	Domain.CheckParams();

	//Domain.CheckIntrExtr(&SCOX);
	//Domain.CheckRetraction(&SCOX);
	//Domain.CheckDiffRetraction(&SCOX);
	//Domain.CheckIsometryofVectorTransport(&SCOX);
	//Domain.CheckIsometryofInvVectorTransport(&SCOX);
	//Domain.CheckLockingCondition(&SCOX);
	//Domain.CheckcoTangentVector(&SCOX);
	//Domain.CheckVecTranComposeInverseVecTran(&SCOX);
	//Domain.CheckTranHInvTran(&SCOX);
	//Domain.CheckHaddScaledRank1OPE(&SCOX);

	NSOLyapunov Prob(A, M, C, n, p);
	/*The domain of the problem is a Stiefel manifold*/
	Prob.SetDomain(&Domain);

	Prob.CheckGradHessian(&SCOX);


	LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, &SCOX);
	LRBFGSsolver->Debug = ITERRESULT; //ITERRESULT;// 
	LRBFGSsolver->Max_Iteration = 4000;
	LRBFGSsolver->Tolerance = 1e-7;
	LRBFGSsolver->OutputGap = 50;
	LRBFGSsolver->LengthSY = 4;
	LRBFGSsolver->CheckParams();
	LRBFGSsolver->Run();

	Prob.CheckGradHessian(LRBFGSsolver->GetXopt());

	delete LRBFGSsolver;

	//RCG *RCGsolver = new RCG(&Prob, &SCOX);
	//RCGsolver->Debug = ITERRESULT; //ITERRESULT;// 
	//RCGsolver->Max_Iteration = 4000;
	//RCGsolver->Tolerance = 1e-7;
	//RCGsolver->OutputGap = 50;
	//RCGsolver->Run();

	//delete RCGsolver;

	delete[] A;
}


#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 7)
	{
		mexErrMsgTxt("The number of arguments should be at least seven.\n");
	}
	double *A, *M, *C, *X;
	A = mxGetPr(prhs[0]);
	M = mxGetPr(prhs[1]);
	C = mxGetPr(prhs[2]);
	X = mxGetPr(prhs[3]);
	/* dimensions of input matrices */
	integer n, p, HasHHR, metric;
	n = mxGetM(prhs[0]);
	p = mxGetN(prhs[3]);
	if (mxGetN(prhs[0]) != n || mxGetM(prhs[1]) != n || mxGetN(prhs[1]) != n || mxGetM(prhs[2]) != n || mxGetN(prhs[2]) != n || mxGetM(prhs[3]) != n)
	{
		mexErrMsgTxt("The size of A or the size of M or the size of C or the size of X is not correct.\n");
	}
	HasHHR = static_cast<integer> (mxGetScalar(prhs[4]));
	metric = static_cast<integer> (mxGetScalar(prhs[5]));

	bool isAsparse = mxIsSparse(prhs[0]);
	integer nzmaxA;
	size_t *irA = nullptr, *jcA = nullptr;
	int *inirA = nullptr, *injcA = nullptr;
	if (isAsparse)
	{
		nzmaxA = mxGetNzmax(prhs[1]);
		irA = mxGetIr(prhs[1]);
		jcA = mxGetJc(prhs[1]);

		inirA = new int[2 * nzmaxA];
		injcA = inirA + nzmaxA;
		for (integer i = 0; i < n; i++)
		{
			for (unsigned long long j = jcA[i]; j < jcA[i + 1]; j++)
			{
				/*row: ir[j], column: i, entry: A[j]*/
				inirA[j] = irA[j];
				injcA[j] = i;
			}
		}
	}
	
	genrandseed(0);
	CheckMemoryDeleted = new std::map<integer *, integer>;

	// Obtain an initial iterate
	NSOVariable NSOX(n, p);
	double *NSOXptr = NSOX.ObtainWriteEntireData();
	for (integer i = 0; i < n * p; i++)
		NSOXptr[i] = X[i];

	NSOVariable *NSOsoln = nullptr;
	if (nrhs >= 8)
	{
		double *soln = mxGetPr(prhs[7]); /*soln: n by p*/
		NSOsoln = new NSOVariable(n, p);
		double *NSOsolnptr = NSOsoln->ObtainWriteEntireData();
		for (integer i = 0; i < n * p; i++)
		{
			NSOsolnptr[i] = soln[i];
		}
	}

	// Define the manifold
	NStQOrth Domain(n, p);

	// Define the matrix completion problem
	NSOLyapunov Prob(A, M, C, n, p);
	Prob.SetDomain(&Domain);

	Domain.SetHasHHR(HasHHR != 0);
	if (metric == 1)
		Domain.ChooseNSOParamsSet1();
	if (metric == 2)
		Domain.ChooseNSOParamsSet2();
	if (metric == 3)
		Domain.ChooseNSOParamsSet3();
	if (metric == 4)
		Domain.ChooseNSOParamsSet4();
	if (metric == 5)
		Domain.ChooseNSOParamsSet5();
	if (metric == 6)
		Domain.ChooseNSOParamsSet6();
	//Domain.CheckParams();

	// Call the function defined in DriverMexProb.h
	ParseSolverParamsAndOptimizing(prhs[6], &Prob, &NSOX, NSOsoln, plhs);

	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	delete NSOsoln;
	delete[] inirA;

	return;
}

#endif