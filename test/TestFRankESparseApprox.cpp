#include "test/TestFRankESparseApprox.h"

using namespace ROPTLIB;

void testFRankESparseApprox(void)
{
	//	FILE *ttt = nullptr;
	//	freopen_s(&ttt, "./log.txt", "w", stdout);
	integer m = 10, n = 10, r = 2;
	//integer m = 100, n = 15, r = 5;
	FixedRankE Domain(m, n, r);
	Domain.SetHasHHR(false);
    Variable InitialX = Domain.RandominManifold();
//    InitialX.Print("X:", false);//---
//	Variable InitialX(m, n, r);
//	InitialX.RandInManifold();
//  Domain.CheckParams();
//  Domain.CheckIntrExtr(&InitialX);
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

	//InitialX.Print("initialX:");

	// Generate the matrices in the Low rank approximation problem.
    Vector A(m, n); A.RandGaussian();
    realdp lambda = 1;
    integer lengthW = 1;
    
    FRankESparseApprox Prob(A, lambda, m, n, r, lengthW);
	Prob.SetDomain(&Domain);

//	Prob.CheckGradHessian(InitialX);
    
    ManPG *ManPGsolver = new ManPG(&Prob, &InitialX);
    ManPGsolver->Max_Iteration = 200;
//    ManPGsolver->Variant = LSPG_REGULAR; //-- LSPG_REGULAR; //-- LSPG_ADALIPSCHITZ;
    ManPGsolver->CheckParams();
    ManPGsolver->Run();
//    ManPGsolver->GetXopt().Print("xopt");//---
    delete ManPGsolver;
};

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{/*input: A, lambda, initX, SolverParams*/
    if(nrhs < 4)
    {
        mexErrMsgTxt("The number of arguments should be at least four.\n");
    }
    realdp *A;
    realdp lambda;
    A = mxGetPr(prhs[0]);
    lambda = static_cast<realdp> (mxGetScalar(prhs[1]));
    /* dimensions of input matrices */
    integer m, n, r, HasHHR;
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);

    Vector InitialX(m, n);
    mexProblem::ObtainElementFromMxArray(&InitialX, prhs[2]);
    r = InitialX.Field("U").Getcol();
    
    // Generate the matrices in the Low rank approximation problem.
    Vector AA(m, n);
    realdp *AAptr = AA.ObtainWriteEntireData();
    for(integer i = 0; i < m * n; i++)
        AAptr[i] = A[i];
    
    HasHHR = 0;

    genrandseed(0);

    CheckMemoryDeleted = new std::map<integer *, integer>;
    
    // Define the manifold
    FixedRankE Domain(m, n, r);
    
    integer lengthW = 1;
    
    FRankESparseApprox Prob(AA, lambda, m, n, r, lengthW);
    
    Prob.SetDomain(&Domain);
    
    Domain.SetHasHHR(HasHHR != 0);

    // Call the function defined in DriverMexProb.h
    ParseSolverParamsAndOptimizing(prhs[3], &Prob, &InitialX, plhs);

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
