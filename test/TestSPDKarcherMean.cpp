
#include "test/TestSPDKarcherMean.h"

using namespace ROPTLIB;

void testSPDKarcherMean(void)
{
	/*Randomly generate a point on the SPD manifold*/
	integer n = 3, num = 4;

	// Define the manifold
	SPDManifold Domain(n);
//	Domain.SetHasHHR(true); /*set whether the manifold uses the idea in [HGA2015, Section 4.3] or not*/
    Domain.ChooseParamsSet1();
    Variable SPDX = Domain.RandominManifold();

    Vector EE(n, n), tmp(n, n);
    Vector Ls(1, &EE, num);
    Ls.NewMemoryOnWrite();
    for(integer i = 0; i < num; i++)
    {
        tmp.RandGaussian(); tmp.QRDecom();
        Ls.GetElement(i) = tmp.Field("_R").GetTranspose();
    }
    
	// Define the problem
	SPDKarcherMean Prob(Ls, n, num);
	/*The domain of the problem is a SPD manifold*/
	Prob.SetDomain(&Domain);
    Domain.CheckParams();
//    Prob.SetNumGradHess(true);
	Prob.CheckGradHessian(SPDX);
//    return;
	/*Output the parameters of the domain manifold*/
	//Domain.CheckParams();

	/*Check the correctness of the manifold operations*/
//	Domain.CheckIntrExtr(SPDX);
//	Domain.CheckRetraction(SPDX);
//	Domain.CheckDiffRetraction(SPDX);
//	Domain.CheckLockingCondition(SPDX);
	//Domain.CheckcoTangentVector(&SPDX);
	//Domain.CheckIsometryofVectorTransport(&SPDX);
	//Domain.CheckIsometryofInvVectorTransport(&SPDX);
	//Domain.CheckVecTranComposeInverseVecTran(&SPDX);
	//Domain.CheckTranHInvTran(&SPDX);
	//Domain.CheckHaddScaledRank1OPE(&SPDX);
    
	// test LRBFGS
	//printf("********************************Test Geometric mean in LRBFGS*************************************\n");
	RSD *LRBFGSsolver = new RSD(&Prob, &SPDX);
	LRBFGSsolver->LineSearch_LS = LSSM_ARMIJO;
	LRBFGSsolver->Verbose = FINALRESULT; //ITERRESULT;//
	LRBFGSsolver->Max_Iteration = 50;
	LRBFGSsolver->Tolerance = static_cast<realdp> (1e-6);
	LRBFGSsolver->Accuracy = static_cast<realdp> (1e-4);
	LRBFGSsolver->Finalstepsize = 1;
    LRBFGSsolver->LineSearch_LS = LSSM_STRONGWOLFE;
	LRBFGSsolver->CheckParams();
	LRBFGSsolver->Run();
	if (LRBFGSsolver->Getnormgfgf0() < 1e-6)
		printf("SUCCESS!\n");
	else
		printf("FAIL!\n");
	//// Check gradient and Hessian
	//Prob.CheckGradHessian(&SPDX);
//	Prob.CheckGradHessian(LRBFGSsolver->GetXopt());

	//LRBFGS *LRBFGSsolver2 = new LRBFGS(&Prob, &SPDX, LRBFGSsolver->GetXopt());
	//LRBFGSsolver2->LineSearch_LS = ARMIJO;
	//LRBFGSsolver2->Debug = ITERRESULT; //ITERRESULT;// 
	//LRBFGSsolver2->Max_Iteration = 500;
	//LRBFGSsolver2->Tolerance = 1e-5;
	//LRBFGSsolver2->Accuracy = 1e-4;
	//LRBFGSsolver2->Finalstepsize = 1;
	//LRBFGSsolver2->CheckParams();
	//LRBFGSsolver2->Run();
	//realdp *dists = LRBFGSsolver2->GetdistSeries();
	//ForDebug::Print("Dist:", dists, LRBFGSsolver2->GetlengthSeries());
	//delete LRBFGSsolver2;
	delete LRBFGSsolver;

};

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 5)
	{
		mexErrMsgTxt("The number of arguments should be at least five.\n");
	}
	realdp *Ls, *X;
	integer n, N, HasHHR, ParamSet;
	Ls = mxGetPr(prhs[0]);
	X = mxGetPr(prhs[1]);
	n = mxGetM(prhs[1]);
	const mwSize *ptrdims = mxGetDimensions(prhs[0]);
	if (mxGetNumberOfDimensions(prhs[0]) == 2)
		N = 1;
	else
		N = ptrdims[2];
    
    if (ptrdims[1] != n || ptrdims[0] != n)
    {
        mexErrMsgTxt("The size of matrix C is not correct.\n");
    }
    if (mxGetM(prhs[1]) != n || mxGetN(prhs[1]) != n)
    {
        mexErrMsgTxt("The size of the initial X is not correct!\n");
    }

    HasHHR = static_cast<integer> (mxGetScalar(prhs[2]));
    ParamSet = static_cast<integer> (mxGetScalar(prhs[3]));

    Vector EE(n, n), tmp(n, n);
    Vector LLs(1, &EE, N);
    realdp *LLsptr = LLs.ObtainWriteEntireData();
    for(integer i = 0; i < n * n * N; i++)
        LLsptr[i] = Ls[i];
    
	genrandseed(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;
	//testStieSoftICA(Cs, n, p, N, X, Xopt);

	// Define the manifold
	SPDManifold Domain(n);

	if (ParamSet == 1)
		Domain.ChooseParamsSet();
	else
    if (ParamSet == 2)
		Domain.ChooseParamsSet2();
    else
    if (ParamSet == 3)
        Domain.ChooseParamsSet3();
    else
        Domain.ChooseParamsSet4();
    
    Variable initX = Domain.RandominManifold();
    realdp *initXptr = initX.ObtainWriteEntireData();
    for(integer i = 0; i < n * n; i++)
        initXptr[i] = X[i];
    
    
    SPDKarcherMean Prob(LLs, n, N);
    Prob.SetDomain(&Domain);
	Domain.SetHasHHR((HasHHR != 0));
    
	ParseSolverParamsAndOptimizing(prhs[4], &Prob, &initX, plhs);
    
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
