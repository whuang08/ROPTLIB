
#include "test/TestSphereSparsestVector.h"

using namespace ROPTLIB;

void testSphereSparsestVector(void)
{
	// size of the matrix Q
	integer m = 100, n =10;

	// Generate the matrix
    Vector Q(m, n);
    Q.RandGaussian();

//    Variable SphereX(n);
//    SphereX.RandInManifold();

    // Define the manifold
    Sphere Domain(n);
    Domain.SetHasHHR(false);
//    Domain.ChooseParamsSet2();
    //Domain.SetHasHHR(true); /*set whether the manifold uses the idea in [HGA2015, Section 4.3] or not*/
    Variable SphereX = Domain.RandominManifold();
    
    // Define the SparestVector problem
    SphereSparsestVector Prob(Q);
    /*The domain of the problem is a Stiefel manifold*/
    Prob.SetDomain(&Domain);

    /*Output the parameters of the domain manifold*/
    Domain.CheckParams();

    //Prob.CheckGradHessian(&SphereX);
    
    //Domain.CheckRetraction(&StieX);
    //Domain.CheckDiffRetraction(&StieX);
    //Domain.CheckLockingCondition(&StieX);
    //Domain.CheckcoTangentVector(&StieX);
    //Domain.CheckIsometryofVectorTransport(&StieX);
    //Domain.CheckIsometryofInvVectorTransport(&StieX);
    //Domain.CheckVecTranComposeInverseVecTran(&StieX);
    //Domain.CheckTranHInvTran(&StieX);
    //Domain.CheckHaddScaledRank1OPE(&StieX);
    
//    LRBFGSSub *LRBFGSSubsolver = new LRBFGSSub(&Prob, &SphereX);
    RGS *LRBFGSSubsolver = new RGS(&Prob, &SphereX);
//    RGS *LRBFGSSubsolver = new RGS(&Prob, &SphereX);
    LRBFGSSubsolver->Verbose = ITERRESULT;//--- FINALRESULT;
    LRBFGSSubsolver->OutputGap = 1;
//    LRBFGSSubsolver->lambdaLower = static_cast<realdp> (1e-3);
//    LRBFGSSubsolver->lambdaUpper = static_cast<realdp> (1e3);
    LRBFGSSubsolver->Tolerance = static_cast<realdp> (1e-6);
    LRBFGSSubsolver->Max_Iteration = 200;
    //LRBFGSSubsolver->Stop_Criterion = FUN_REL;
    LRBFGSSubsolver->CheckParams();
    LRBFGSSubsolver->Run();
    if (LRBFGSSubsolver->Getnormndnd0() < 1e-6)
        printf("SUCCESS!\n");
    else
        printf("FAIL!\n");
    delete LRBFGSSubsolver;

    //RBFGSSubsolver = new RBFGSSub(&Prob, &StieX);
    //RBFGSSubsolver->Debug = FINALRESULT;
    //RBFGSSubsolver->OutputGap = 10;
    //RBFGSSubsolver->lambdaLower = 1e-7;
    //RBFGSSubsolver->lambdaUpper = 1e7;
    //RBFGSSubsolver->CheckParams();
    //RBFGSSubsolver->Run();
    //delete RBFGSSubsolver;
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
	if(nrhs < 5)
	{
		mexErrMsgTxt("The number of arguments should be at least five.\n");
	}
	realdp *Q, *X;
	Q = mxGetPr(prhs[0]);
	X = mxGetPr(prhs[1]);
	/* dimensions of input matrices */
	integer m, n, HasHHR, Paramset;
	m = mxGetM(prhs[0]);
	n = mxGetM(prhs[1]);

	/*Check the correctness of the inputs*/
	if(mxGetN(prhs[0]) != n)
	{
		mexErrMsgTxt("The size of matrix Q or vector X is not correct.\n");
	}
	HasHHR = static_cast<integer> (mxGetScalar(prhs[2]));
	Paramset = static_cast<integer> (mxGetScalar(prhs[3]));

	genrandseed(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;

    Vector QQ(m, n);
    realdp *QQptr = QQ.ObtainWriteEntireData();
    for(integer i = 0; i < m * n; i++)
        QQptr[i] = Q[i];
    Sphere Domain(n);
    Variable initX = Domain.RandominManifold();
    realdp *initXptr = initX.ObtainWriteEntireData();
    for(integer i = 0; i < n; i++)
        initXptr[i] = X[i];
    
    Domain.SetHasHHR(HasHHR);
    
	if (Paramset == 1)
		Domain.ChooseParamsSet1();
	else if (Paramset == 2)
		Domain.ChooseParamsSet2();
	else if (Paramset == 3)
		Domain.ChooseParamsSet3();
	else if (Paramset == 4)
		Domain.ChooseParamsSet4();

	// Define the SparestVector problem
	SphereSparsestVector Prob(QQ);
	Prob.SetDomain(&Domain);

	Domain.SetHasHHR(HasHHR != 0);
	//Domain.CheckParams();

	// Call the function defined in DriverMexProb.h
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
