#include "test/TestFRankEWeightApprox.h"

using namespace ROPTLIB;

void testFRankEWeightApprox(void)
{
	//	FILE *ttt = nullptr;
	//	freopen_s(&ttt, "./log.txt", "w", stdout);
	integer m = 5, n = 4, r = 2;
	//integer m = 100, n = 15, r = 5;
	FixedRankE Domain(m, n, r);
	Domain.SetHasHHR(true);
    Variable InitialX = Domain.RandominManifold();
//    InitialX.Print("X:", false);//---
//	Variable InitialX(m, n, r);
//	InitialX.RandInManifold();
	//Domain.CheckParams();
	//Domain.CheckIntrExtr(&InitialX);
//	Domain.CheckRetraction(InitialX);
//    Domain.CheckDiffRetraction(InitialX, true);
//    Domain.CheckLockingCondition(InitialX);
//	Domain.CheckcoTangentVector(&InitialX);
//	Domain.CheckDiffRetraction(&InitialX, false);
//	Domain.CheckIsometryofVectorTransport(&InitialX);
//
//	Domain.CheckLockingCondition(&InitialX);
//	Domain.CheckIsometryofInvVectorTransport(&InitialX);
//	Domain.CheckVecTranComposeInverseVecTran(&InitialX);
//	Domain.CheckTranHInvTran(&InitialX);

//	InitialX.Print("initialX:");

	// Generate the matrices in the Low rank approximation problem.
    Vector A(m, n); A.RandGaussian();
    Vector O(m * n, m * n);
    O.RandGaussian();
    O = O.GetOrth();
    Vector D(m * n);
    D.RandUnform();
    D = D + 0.1;
    Vector W = O.GetTranspose() * D.GetDiagTimesM(O);
    
	FRankEWeightApprox Prob(A, W, m, n, r);
	Prob.SetDomain(&Domain);

	Prob.CheckGradHessian(InitialX);

	LRBFGS *RSDsolver = new LRBFGS(&Prob, &InitialX);
	//->LineSearch_LS = ARMIJO;
	//RSDsolver->LS_beta = 0.01;
	//RSDsolver->RCGmethod = DAI_YUAN;
	RSDsolver->Verbose = FINALRESULT;
	RSDsolver->OutputGap = 1;
	RSDsolver->Max_Iteration = 70;
	RSDsolver->InitSteptype = LSSM_ONESTEP;
	//RSDsolver->CheckParams();
	//RSDsolver->Accuracy = 1e-6;
	RSDsolver->Tolerance = static_cast<realdp> (1e-6);
	RSDsolver->Run();
	if (RSDsolver->Getnormgfgf0() < 1e-6)
		printf("SUCCESS!\n");
	else
		printf("FAIL!\n");
	//Prob.CheckGradHessian(&InitialX);//--
	//Prob.CheckGradHessian(RSDsolver->GetXopt());//--
	//RSDsolver->GetXopt()->Print("Xopt");
    
    Prob.CheckGradHessian(RSDsolver->GetXopt());
    
    Vector xopt = RSDsolver->GetXopt();
    Prob.MinMaxEigValHess(InitialX).Print("HessEigVal:");
    Prob.MinMaxEigValHess(xopt).Print("HessEigVal:");
    
	delete RSDsolver;
};

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs < 5)
    {
        mexErrMsgTxt("The number of arguments should be at least five.\n");
    }
    realdp *A, *W;
    A = mxGetPr(prhs[0]);
    W = mxGetPr(prhs[1]);
    /* dimensions of input matrices */
    integer m, n, r, HasHHR;
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);

    /*Check the correctness of the inputs*/
    if(mxGetN(prhs[1]) != m * n || mxGetM(prhs[1]) != m * n)
    {
        mexErrMsgTxt("The size of matrix W is not correct.\n");
    }
    if (!mxIsStruct(prhs[2]))
    {
        mexErrMsgTxt("The third argument, initial iterate, is not a structure.");
    }
    Vector InitialX(m, n);
    mexProblem::ObtainElementFromMxArray(&InitialX, prhs[2]);
    r = InitialX.Field("U").Getcol();
    
    HasHHR = static_cast<integer> (mxGetScalar(prhs[3]));

    genrandseed(0);

    CheckMemoryDeleted = new std::map<integer *, integer>;
    
    // Define the manifold
    FixedRankE Domain(m, n, r);
    
    
    // Generate the matrices in the Low rank approximation problem.
    Vector AA(m, n);
    realdp *AAptr = AA.ObtainWriteEntireData();
    for(integer i = 0; i < m * n; i++)
        AAptr[i] = A[i];
    
    Vector WW(m * n, m * n);
    realdp *WWptr = WW.ObtainWriteEntireData();
    for(integer i = 0; i < m * n * m * n; i++)
        WWptr[i] = W[i];
    
    FRankEWeightApprox Prob(AA, WW, m, n, r);
    Prob.SetDomain(&Domain);
    
    Domain.SetHasHHR(HasHHR != 0);

    // Call the function defined in DriverMexProb.h
    ParseSolverParamsAndOptimizing(prhs[4], &Prob, &InitialX, plhs);

    std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
    for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
    {
        if (iter->second != 1)
            printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
    }
    delete CheckMemoryDeleted;
    return;
}

#endif
