
#include "test/TestGrassRQ.h"

using namespace ROPTLIB;

void testGrassRQ(void)
{
    unsigned tt = (unsigned)time(NULL);
    tt = 2; /*The following test is only for random seed zero*/
    std::cout << "seed SB:" << tt << std::endl;//---
    genrandseed(tt);

    // size of the Stiefel manifold
    integer n = 5, p = 2;
    // Generate the matrices in the Brockett problem.
    Vector B(n, n);
    B.RandGaussian();
    B = B + B.GetTranspose();

    // Define the manifold
    Grassmann Domain(n, p);
    Variable GrassX = Domain.RandominManifold();

                            // Define the Brockett problem
    GrassRQ Prob(B, n, p);
    /*The domain of the problem is a Stiefel manifold*/
    Prob.SetDomain(&Domain);
    Domain.CheckParams();

//    Domain.CheckIntrExtr(GrassX);
//    Domain.CheckRetraction(GrassX);
//    Domain.CheckDiffRetraction(GrassX, true);
//    Domain.CheckLockingCondition(StieX);
//    Domain.CheckcoTangentVector(StieX);
//    Domain.CheckIsometryofVectorTransport(StieX);
//    Domain.CheckIsometryofInvVectorTransport(StieX);
//    Domain.CheckVecTranComposeInverseVecTran(StieX);
//    Domain.CheckTranHInvTran(StieX);
//    Domain.CheckHaddScaledRank1OPE(StieX);
//    return;
    Prob.CheckGradHessian(GrassX);

    LRTRSR1 *RSDsolver = new LRTRSR1(&Prob, &GrassX);
    RSDsolver->Verbose = FINALRESULT;//--- FINALRESULT;
    RSDsolver->Max_Iteration = 100;
    RSDsolver->OutputGap = 1;
    RSDsolver->CheckParams();
    RSDsolver->Run();

    Prob.CheckGradHessian(RSDsolver->GetXopt());

    Vector xopt = RSDsolver->GetXopt();
    Prob.MinMaxEigValHess(GrassX).Print("EigValHess:");
    Prob.MinMaxEigValHess(xopt).Print("EigValHess:");
    delete RSDsolver;
}

/*If it is compiled in Matlab, then the following "mexFunction" is used as the entrance.*/
#ifdef MATLAB_MEX_FILE

/*Help to check the memory leakage problem. No necesary any more.*/
std::map<integer *, integer> *CheckMemoryDeleted;

/*This function checks the number and formats of input parameters.
nlhs: the number of output in mxArray format
plhs: the output objects in mxArray format
nrhs: the number of input in mxArray format
prhs: the input objects in mxArray format */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if(nrhs < 4)
	{
		mexErrMsgTxt("The number of arguments should be at least four.\n");
	}
	realdp *B, *X, *Xopt, *soln;
	B = mxGetPr(prhs[0]);
	X = mxGetPr(prhs[1]);
	/* dimensions of input matrices */
	integer p, n, HasHHR;
	n = mxGetM(prhs[0]);
	p = mxGetN(prhs[1]);

	/*Check the correctness of the inputs*/
	if(mxGetN(prhs[0]) != n)
	{
		mexErrMsgTxt("The size of matrix B is not correct.\n");
	}
	if(mxGetM(prhs[1]) != n)
	{
		mexErrMsgTxt("The size of the initial X is not correct!\n");
	}
	HasHHR = static_cast<integer> (mxGetScalar(prhs[2]));

	printf("(n, p):%d,%d\n", n, p);

	/*create output matrix*/
	plhs[0] = mxCreateDoubleMatrix(n, p, mxREAL);
	Xopt = mxGetPr(plhs[0]);

	genrandseed(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;

	// Define the manifold
	Grassmann Domain(n, p);
    
    Variable initX = Domain.RandominManifold();
    realdp *initXptr = initX.ObtainWriteEntireData();
    for(integer i = 0; i < n * p; i++)
        initXptr[i] = X[i];

	// Define the Rayleigh Quotient problem
    Vector BB(n, n);
    realdp *BBptr = BB.ObtainWriteEntireData();
    for(integer i = 0; i < n * n; i++)
        BBptr[i] = B[i];
	GrassRQ Prob(BB, n, p);
	Prob.SetDomain(&Domain);

	Domain.SetHasHHR(HasHHR != 0);

	/* Call the function defined in DriverMexProb.h */
	ParseSolverParamsAndOptimizing(prhs[3], &Prob, &initX, plhs);

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
