#include "test/TestSFRQLyapunov.h"

using namespace ROPTLIB;

void testSFRQLyapunov(void)
{
	integer n = 10, p = 2, pC = 1;

	integer nzmaxA = n + 2 * (n - 1);
	realdp *A = new realdp[n + 2 * (n - 1)];
	integer *inirA = new integer[2 * nzmaxA];
	integer *injcA = inirA + nzmaxA;
	for (integer i = 0; i < n; i++)
	{
		A[i] = 2;
		inirA[i] = i;
		injcA[i] = i;
	}
	for (integer i = 0; i < n - 1; i++)
	{
		A[n + i] = -1;
		inirA[n + i] = i;
		injcA[n + i] = i + 1;
		A[2 * n - 1 + i] = -1;
		inirA[2 * n - 1 + i] = i + 1;
		injcA[2 * n - 1 + i] = i;
	}
    SparseMatrix sA(n, n, inirA, injcA, A, nzmaxA);

    integer nzmaxM = n;
    realdp *M = new realdp[n];
    integer *inirM = new integer[2 * nzmaxM];
    integer *injcM = inirM + nzmaxM;
    for(integer i = 0; i < n; i++)
    {
        M[i] = 1;
        inirM[i] = i;
        injcM[i] = i;
    }
    SparseMatrix sM(n, n, inirM, injcM, M, nzmaxM);
    
    Vector C(n, pC);
    C.RandGaussian();
    
    SymFixedRankQ Domain(n, p);
    Vector InitX = Domain.RandominManifold();
    Domain.ChooseParamsSet2();
    Domain.CheckParams();
    
//    Domain.CheckIntrExtr(InitX);
//    Domain.CheckRetraction(InitX);
//    Domain.CheckcoTangentVector(InitX);
//    Domain.CheckDiffRetraction(InitX, false);
//    Domain.CheckIsometryofVectorTransport(InitX);
//
//    Domain.CheckLockingCondition(InitX);
//    Domain.CheckIsometryofInvVectorTransport(InitX);
//    Domain.CheckVecTranComposeInverseVecTran(InitX);
//    Domain.CheckTranHInvTran(InitX);
//    return;
    
    SFRQLyapunov Prob(sA, sM, C, p);
    
    Prob.SetDomain(&Domain);

    Prob.CheckGradHessian(InitX);

    LRBFGS *RSDsolver = new LRBFGS(&Prob, &InitX);
//    RSD *RSDsolver = new RSD(&Prob, &InitX);
    //->LineSearch_LS = ARMIJO;
    //RSDsolver->LS_beta = 0.01;
    //RSDsolver->RCGmethod = DAI_YUAN;
    RSDsolver->Verbose = ITERRESULT;//-- FINALRESULT;
    RSDsolver->OutputGap = 1;
    RSDsolver->Max_Iteration = 500;
    RSDsolver->Accuracy = 1e-6;
    RSDsolver->Finalstepsize = 1;
    RSDsolver->Tolerance = 1e-6;
//    RSDsolver->LengthSY = 0;
//    RSDsolver->nu = 0;
    RSDsolver->LS_ratio1 = 0.3;
    RSDsolver->LS_ratio2 = 0.3;
    RSDsolver->InitSteptype = LSSM_ONESTEP;
    RSDsolver->CheckParams();
    RSDsolver->Run();
    //Prob.CheckGradHessian(&InitialX);//--
    //Prob.CheckGradHessian(RSDsolver->GetXopt());//--
    if (RSDsolver->Getnormgfgf0() < 1e-6)
        printf("SUCCESS!\n");
    else
        printf("FAIL!\n");

    Prob.CheckGradHessian(RSDsolver->GetXopt());

    Prob.MinMaxEigValHess(RSDsolver->GetXopt()).Print("eigs:");//--
    delete RSDsolver;
                      
    delete [] A;
    delete [] M;
    delete [] inirA;
    delete [] inirM;
}


#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 7)
	{
		mexErrMsgTxt("The number of arguments should be at least seven.\n");
	}
	realdp *A, *M, *C, *X;
	A = mxGetPr(prhs[0]);
	M = mxGetPr(prhs[1]);
	C = mxGetPr(prhs[2]);
	X = mxGetPr(prhs[3]);
	/* dimensions of input matrices */
	integer n, p, pC, HasHHR, ParamSet;
	n = mxGetM(prhs[0]);
	p = mxGetN(prhs[3]);
	pC = mxGetN(prhs[2]);
	if (mxGetN(prhs[0]) != n || mxGetM(prhs[1]) != n || mxGetN(prhs[1]) != n || mxGetM(prhs[2]) != n || mxGetM(prhs[3]) != n)
	{
		mexErrMsgTxt("The size of A or the size of M or the size of C or the size of X is not correct.\n");
	}
	HasHHR = static_cast<integer> (mxGetScalar(prhs[4]));
	ParamSet = static_cast<integer> (mxGetScalar(prhs[5]));

	bool isAsparse = mxIsSparse(prhs[0]);
	integer nzmaxA;
	size_t *irA = nullptr, *jcA = nullptr;
	integer *inirA = nullptr, *injcA = nullptr;
	if (! isAsparse)
	{
        mexErrMsgTxt("A need be a sparse matrix.\n");
	}
    nzmaxA = mxGetNzmax(prhs[0]);
    irA = mxGetIr(prhs[0]);
    jcA = mxGetJc(prhs[0]);

    inirA = new integer[2 * nzmaxA];
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
    
    SparseMatrix sA(n, n, inirA, injcA, A, nzmaxA);
    
	bool isMsparse = mxIsSparse(prhs[1]);
	integer nzmaxM = 0;
	size_t *irM = nullptr, *jcM = nullptr;
	integer *inirM = nullptr, *injcM = nullptr;
	if (! isMsparse)
	{
        mexErrMsgTxt("M need be a sparse matrix.\n");
	}
    nzmaxM = mxGetNzmax(prhs[1]);
    irM = mxGetIr(prhs[1]);
    jcM = mxGetJc(prhs[1]);

    inirM = new integer[2 * nzmaxM];
    injcM = inirM + nzmaxM;
    for (integer i = 0; i < n; i++)
    {
        for (unsigned long long j = jcM[i]; j < jcM[i + 1]; j++)
        {
            /*row: ir[j], column: i, entry: M[j]*/
            inirM[j] = irM[j];
            injcM[j] = i;
        }
    }
    
    SparseMatrix sM(n, n, inirM, injcM, M, nzmaxM);
    
	genrandseed(0);
	CheckMemoryDeleted = new std::map<integer *, integer>;

	// Define the manifold
	SymFixedRankQ Domain(n, p);
    Vector initialX = Domain.RandominManifold();
    realdp *initialXptr = initialX.ObtainWriteEntireData();
    for(integer i = 0; i < n * p; i++)
        initialXptr[i] = X[i];
    
    Vector CC(n, pC);
    realdp *CCptr = CC.ObtainWriteEntireData();
    for(integer i = 0; i < n * pC; i++)
        CCptr[i] = C[i];
	/* Define the matrix completion problem */
    SFRQLyapunov Prob(sA, sM, CC, p);
	Prob.SetDomain(&Domain);

	Domain.SetHasHHR(HasHHR != 0);
	if (ParamSet == 1)
		Domain.ChooseParamsSet1();
	if (ParamSet == 2)
		Domain.ChooseParamsSet2();
	if (ParamSet == 3)
		Domain.ChooseParamsSet3();
	if (ParamSet == 4)
		Domain.ChooseParamsSet4();
	if (ParamSet == 5)
		Domain.ChooseParamsSet5();
	if (ParamSet == 6)
		Domain.ChooseParamsSet6();

	/* Call the function defined in DriverMexProb.h */
	ParseSolverParamsAndOptimizing(prhs[6], &Prob, &initialX, plhs);

	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	delete[] inirA;
	delete[] inirM;

	return;
}

#endif
