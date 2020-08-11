#include "test/TestCSFRQPhaseRetrieval.h"

#ifdef ROPTLIB_WITH_FFTW

using namespace ROPTLIB;

void testCSFRQPhaseRetrieval(void)
{
    integer n1 = 128, n2 = 128, r = 2, l = 6;
	integer n = n1 * n2, m = n * l;
    realdp kappa = 1e-8;

	CSymFixedRankQ Domain(n, r);
    Vector InitialX = Domain.RandominManifold();
    Domain.ChooseParamsSet1();
//	Domain.CheckParams();
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

	// Generate the matrices in the Low rank approximation problem.
//    Vector b(m); b.RandGaussian();
    Vector masks(n, l, "complex"); masks.RandGaussian();
    
    Vector xtrue(n, r, "complex"); xtrue.RandGaussian();

    Vector ZY(m, r, "complex");
    realdp sqn = sqrt(static_cast<realdp> (n1 * n2));
    Vector tmpv;

    for(integer i = 0; i < l; i++)
    {
        for(integer j = 0; j < r; j++)
        {
            tmpv = (masks.GetSubmatrix(0, n - 1, i, i).GetHadamardProduct(xtrue.GetSubmatrix(0, n - 1, j, j))) / sqn;
            ZY.SubmatrixAssignment(i * n, (i + 1) * n - 1, j, j, tmpv.Reshape(n1, n2).GetFFT2D(FFTW_FORWARD).Reshape(n, 1));
        }
    }
    Vector b = ZY.GetTranspose().GetColNormsSquare().GetTranspose();
    
    
    CSFRQPhaseRetrieval Prob(b, masks, kappa, n1, n2, l, r);
	Prob.SetDomain(&Domain);

    Domain.CheckParams();
//	std::cout << "f:" << Prob.f(InitialX) << std::endl;//----
//    Prob.EucGrad(InitialX).Print("EucGrad:");//---

	//Vector *gf = Domain.GetEMPTYINTR()->ConstructEmpty();
	//Prob.Grad(&InitialX, gf);
	//delete gf;
//    Prob.SetNumGradHess(true);

//	Prob.CheckGradHessian(InitialX);
//    return;
//    RTRNewton *RSDsolver = new RTRNewton(&Prob, &InitialX);
    LRBFGS *RSDsolver = new LRBFGS(&Prob, &InitialX);
	//->LineSearch_LS = ARMIJO;
	//RSDsolver->LS_beta = 0.01;
	//RSDsolver->RCGmethod = DAI_YUAN;
    RSDsolver->Verbose = FINALRESULT;
	RSDsolver->OutputGap = 1;
	RSDsolver->Max_Iteration = 100;
	//RSDsolver->CheckParams();
	RSDsolver->Accuracy = 1e-6;
//	RSDsolver->Finalstepsize = 1;
	RSDsolver->Tolerance = 1e-6;
    RSDsolver->CheckParams();
	RSDsolver->Run();
	if (RSDsolver->Getnormgfgf0() < 1e-6)
		printf("SUCCESS!\n");
	else
		printf("FAIL!\n");
//    RSDsolver->GetXopt().Print("xopt:");//---
//	Prob.CheckGradHessian(RSDsolver->GetXopt());//--

	delete RSDsolver;
};

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 6)
	{
		mexErrMsgTxt("The number of arguments should be at least six.\n");
	}
//    std::cout << "h1" << std::endl;//--
	realdp *b, *masks, *X;
	realdp kappa;
    integer HasHHR, ParamSet;
	b = mxGetPr(prhs[0]);
//    std::cout << "h2" << std::endl;//--
    masks = (realdp *) mxGetComplexDoubles(prhs[1]);
//    std::cout << "h3" << std::endl;//--
	X = (realdp *) mxGetComplexDoubles(prhs[2]);
//    std::cout << "h4" << std::endl;//--
    kappa = mxGetScalar(prhs[3]);
//    std::cout << "h5" << std::endl;//--
    ParamSet = mxGetScalar(prhs[4]);
//    std::cout << "h6" << std::endl;//--
    HasHHR = mxGetScalar(prhs[5]);
//    std::cout << "h7" << std::endl;//--
	
	/* dimensions of input matrices */
	integer m, n1, n2, l, n, r;
	const size_t *size = mxGetDimensions(prhs[1]); /* masks */
	n1 = size[0];
	n2 = size[1];
	l = size[2];
	m = mxGetM(prhs[0]);
	n = mxGetM(prhs[2]);
	r = mxGetN(prhs[2]);
	if (n != n1 * n2)
	{
		mexErrMsgTxt("The size of masks or the size of initX is not correct.\n");
	}
	if (m != n * l)
	{
		mexErrMsgTxt("The size of b or the size of masks is not correct.\n");
	}

	genrandseed(0);
	CheckMemoryDeleted = new std::map<integer *, integer>;

    CSymFixedRankQ Domain(n, r);
    if (ParamSet == 1)
        Domain.ChooseParamsSet1();
    else if (ParamSet == 2)
        Domain.ChooseParamsSet2();
    else if (ParamSet == 3)
        Domain.ChooseParamsSet3();
    else
        Domain.ChooseParamsSet4();
    
    Domain.SetHasHHR(HasHHR != 0);
    
    Variable initX = Domain.RandominManifold();
    realdp *initXptr = initX.ObtainWriteEntireData();
    for(integer i = 0; i < initX.Getlength(); i++)
        initXptr[i] = X[i];
    
    Vector mmasks(n, l, "complex");
    realdp *mmasksptr = mmasks.ObtainWriteEntireData();
    for(integer i = 0; i < mmasks.Getlength(); i++)
        mmasksptr[i] = masks[i];
    
    Vector bb(m);
    realdp *bbptr = bb.ObtainWriteEntireData();
    for(integer i = 0; i < bb.Getlength(); i++)
        bbptr[i] = b[i];
    
    CSFRQPhaseRetrieval Prob(bb, mmasks, kappa, n1, n2, l, r);
    Prob.SetDomain(&Domain);
    Prob.SetNumGradHess(false);
    
    // Call the function defined in DriverMexProb.h
    ParseSolverParamsAndOptimizing(prhs[6], &Prob, &initX, plhs);
    
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
#endif
