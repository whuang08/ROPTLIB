#include "test/TestCFRankQ2FBlindDecon2D.h"

#ifdef ROPTLIB_WITH_FFTW

using namespace ROPTLIB;

void testCFRankQ2FBlindDecon2D(void)
{
	integer n1 = 16, n2 = 16, r = 1, L = n1 * n2;

    CFixedRankQ2F Domain(L, L, r);
	Domain.SetHasHHR(false);
	Variable InitialX = Domain.RandominManifold();
	//Domain.CheckParams();

//	Domain.CheckIntrExtr(InitialX);
	//Domain.CheckRetraction(&InitialX);
	//Domain.CheckcoTangentVector(&InitialX);
	//Domain.CheckDiffRetraction(&InitialX, false);
	//Domain.CheckIsometryofVectorTransport(&InitialX);

	//Domain.CheckLockingCondition(&InitialX);
	//Domain.CheckIsometryofInvVectorTransport(&InitialX);
	//Domain.CheckVecTranComposeInverseVecTran(&InitialX);
	//Domain.CheckTranHInvTran(&InitialX);
//	return;

	//InitialX.Print("initialX:");

    Vector y(L, "complex");
    y.RandGaussian();
	// Generate the matrices in the Low rank approximation problem.
    realdp *B = new realdp[L * L * 2 + L * L * 2];
	realdp *C = B + L * L * 2;
	for (integer i = 0; i < L * L * 2 + L * L * 2; i++)
		B[i] = genrandnormal();
	integer nzmaxB = L * L;
	integer nzmaxC = L * L;
	integer *irB = new integer[2 * L * L + 2 * L * L];
	integer *jcB = irB + L * L;
	integer *irC = jcB + L * L;
	integer *jcC = irC + L * L;
	for (integer i = 0; i < L; i++)
	{
		for (integer j = 0; j < L; j++)
		{
			irB[j + i * L] = j;
			jcB[j + i * L] = i;
		}
	}
	for (integer i = 0; i < L; i++)
	{
		for (integer j = 0; j < L; j++)
		{
			irC[j + i * L] = j;
			jcC[j + i * L] = i;
		}
	}
    
    SparseMatrix sB(L, L, irB, jcB, (realdpcomplex *) B, nzmaxB);
    SparseMatrix sC(L, L, irC, jcC, (realdpcomplex *) C, nzmaxC);

    realdp rho = 0;
    CFRankQ2FBlindDecon2D Prob(y, sB, sC, n1, n2, r, rho, 1, 1);
    
	Prob.SetDomain(&Domain);

//	Prob.CheckGradHessian(InitialX);

	LRBFGS *RSDsolver = new LRBFGS(&Prob, &InitialX);
//	RSD *RSDsolver = new RSD(&Prob, &InitialX);
	//->LineSearch_LS = ARMIJO;
	//RSDsolver->LS_beta = 0.01;
	//RSDsolver->RCGmethod = DAI_YUAN;
	RSDsolver->Verbose = FINALRESULT;
	RSDsolver->OutputGap = 1;
	RSDsolver->Max_Iteration = 100;
	RSDsolver->Accuracy = 1e-6;
	RSDsolver->Finalstepsize = 1;
	RSDsolver->Tolerance = 1e-6;
//	RSDsolver->LengthSY = 0;
//	RSDsolver->nu = 0;
	RSDsolver->LS_ratio1 = 0.3;
	RSDsolver->LS_ratio2 = 0.3;
//	RSDsolver->InitSteptype = LSSM_ONESTEP;
	//RSDsolver->CheckParams();
	RSDsolver->Run();
	//Prob.CheckGradHessian(&InitialX);//--
	//Prob.CheckGradHessian(RSDsolver->GetXopt());//--
	if (RSDsolver->Getnormgfgf0() < 1e-6)
		printf("SUCCESS!\n");
	else
		printf("FAIL!\n");
    
//    Prob.CheckGradHessian(RSDsolver->GetXopt());
//    Prob.MinMaxEigValHess(RSDsolver->GetXopt()).Print("eigs:");//--
	delete RSDsolver;

	delete[] B;
	delete[] irB;
};

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ /*TestCFRankQ2FBlindDecon2D(y, B, C, Xinitial, n1, n2, r, HasHHR, SolverParams, rho, d, mu);*/
	if (nrhs < 12)
	{
		mexErrMsgTxt("The number of arguments should be at least twelve.\n");
	}
	realdp *y, *B, *C, *X;
    y = (realdp *) mxGetComplexDoubles(prhs[0]);
    B = (realdp *) mxGetComplexDoubles(prhs[1]);
    C = (realdp *) mxGetComplexDoubles(prhs[2]);

	X = mxGetPr(prhs[3]);
	/* dimensions of input matrices */
	integer L, K, N, HasHHR, n1, n2, r;
	L = mxGetM(prhs[0]);
	K = mxGetN(prhs[1]);
	N = mxGetN(prhs[2]);
	if (K != L || mxGetM(prhs[1]) != L || mxGetM(prhs[2]) != L)
	{
		mexErrMsgTxt("The size of B or the size of C is not correct.\n");
	}
	n1 = static_cast<integer> (mxGetScalar(prhs[4]));
	n2 = static_cast<integer> (mxGetScalar(prhs[5]));
	r = static_cast<integer> (mxGetScalar(prhs[6]));
    HasHHR = static_cast<integer> (mxGetScalar(prhs[7]));
	realdp rho, d, mu;
	rho = mxGetScalar(prhs[9]);
	d = mxGetScalar(prhs[10]);
	mu = mxGetScalar(prhs[11]);
	if (mxGetM(prhs[3]) != 4 * L * r)
	{
		mexErrMsgTxt("The size of initial x is not correct.\n");
	}

	bool isBsparse = mxIsSparse(prhs[1]);
	bool isCsparse = mxIsSparse(prhs[2]);
	integer nzmaxB = 0, nzmaxC = 0;
	size_t *irB = nullptr, *jcB = nullptr, *irC = nullptr, *jcC = nullptr;
	integer *inirB = nullptr, *injcB = nullptr, *inirC = nullptr, *injcC = nullptr;
    
	if (isBsparse)
	{
		nzmaxB = mxGetNzmax(prhs[1]);
		irB = mxGetIr(prhs[1]);
		jcB = mxGetJc(prhs[1]);

		inirB = new integer[2 * nzmaxB];
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

		inirC = new integer[2 * nzmaxC];
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

	// Define the manifold
	CFixedRankQ2F Domain(K, N, r);
    
    SparseMatrix sB(L, K, inirB, injcB, (realdpcomplex *) B, nzmaxB);
    SparseMatrix sC(L, N, inirC, injcC, (realdpcomplex *) C, nzmaxC);
    Vector yy(L, "complex");
    realdp *yyptr = yy.ObtainWriteEntireData();
    for(integer i = 0; i < 2 * L; i++)
        yyptr[i] = y[i];

    CFRankQ2FBlindDecon2D Prob(yy, sB, sC, n1, n2, r, rho, d, mu);
    Prob.SetDomain(&Domain);
    Variable LRX = Domain.RandominManifold();
    realdp *LRXptr = LRX.ObtainWriteEntireData();
    for(integer i = 0; i < LRX.Getlength(); i++)
        LRXptr[i] = X[i];
    
	Domain.SetHasHHR(HasHHR != 0);
//    yy.Print("y:");//--
//    sB.Print("sB:");//---
//    sC.Print("sC:");//---
//    LRX.Print("LRX:");//---
	// Call the function defined in DriverMexProb.h
	ParseSolverParamsAndOptimizing(prhs[8], &Prob, &LRX, plhs);

	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	delete[] inirB;
	delete[] inirC;
//	delete[] y;
//	delete[] B;
//	delete[] C;

	return;
}

#endif
#endif
