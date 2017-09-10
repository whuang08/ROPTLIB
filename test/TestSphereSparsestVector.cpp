
#include "test/TestSphereSparsestVector.h"

using namespace ROPTLIB;

/*If the file is not compiled in Matlab and TESTSPHERESPARESTVECTOR is defined in def.h file, then using the following
main() function as the entrance. */
#if !defined(MATLAB_MEX_FILE) && defined(TESTSPHERESPARSESTVECTOR)

/*Help to check the memory leakage problem. No necesary any more.*/
std::map<integer *, integer> *CheckMemoryDeleted;

int main(void)
{
	/*Set the random seed*/
	unsigned tt = (unsigned)time(NULL);
	tt = 0;
	genrandseed(tt);

	// size of the matrix Q
	integer m = 4, n = 2;

	// Generate the matrix
	double *Q = new double[m * n];
	/*Q is an m by n matrix*/
	for (integer i = 0; i < m * n; i++)
	{
		Q[i] = genrandnormal();
	}

	CheckMemoryDeleted = new std::map<integer *, integer>;
	testSphereSparsestVector(Q, m, n);
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	delete[] Q;
#ifdef _WIN64
#ifdef _DEBUG
	_CrtDumpMemoryLeaks();
#endif
#endif
	return 0;
}

void testSphereSparsestVector(double *Q, integer m, integer n)
{
	SphereVariable SphereX(n);
	SphereX.RandInManifold();

	// Define the manifold
	Sphere Domain(n);
	//Domain.SetHasHHR(true); /*set whether the manifold uses the idea in [HGA2015, Section 4.3] or not*/

	// Define the SparestVector problem
	SphereSparsestVector Prob(Q, m, n);
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

	RBFGSLPSub *RBFGSLPSubsolver = new RBFGSLPSub(&Prob, &SphereX);
	RBFGSLPSubsolver->Debug = ITERRESULT;
	RBFGSLPSubsolver->OutputGap = 1;
	RBFGSLPSubsolver->lambdaLower = 1e-3;
	RBFGSLPSubsolver->lambdaUpper = 1e3;
	RBFGSLPSubsolver->Tolerance = 0.001;
	RBFGSLPSubsolver->Stop_Criterion = FUN_REL;
	RBFGSLPSubsolver->CheckParams();
	RBFGSLPSubsolver->Run();
	delete RBFGSLPSubsolver;
	return;

	//RBFGSLPSubsolver = new RBFGSLPSub(&Prob, &StieX);
	//RBFGSLPSubsolver->Debug = FINALRESULT;
	//RBFGSLPSubsolver->OutputGap = 10;
	//RBFGSLPSubsolver->lambdaLower = 1e-7;
	//RBFGSLPSubsolver->lambdaUpper = 1e7;
	//RBFGSLPSubsolver->CheckParams();
	//RBFGSLPSubsolver->Run();
	//delete RBFGSLPSubsolver;
}

#endif

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
	double *Q, *X, *Xopt;
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
	//	testStieBrockett(B, D, n, p, X, Xopt);

	SphereVariable SphereX(n);
	double *SphereXptr = SphereX.ObtainWriteEntireData();
	for (integer i = 0; i < n; i++)
		SphereXptr[i] = X[i];

	// Define the manifold
	Sphere Domain(n);
	if (Paramset == 1)
		Domain.ChooseSphereParamsSet1();
	else if (Paramset == 2)
		Domain.ChooseSphereParamsSet2();
	else if (Paramset == 3)
		Domain.ChooseSphereParamsSet3();
	else if (Paramset == 4)
		Domain.ChooseSphereParamsSet4();
	else if (Paramset == 5)
		Domain.ChooseSphereParamsSet5();

	// Define the SparestVector problem
	SphereSparsestVector Prob(Q, m, n);
	Prob.SetDomain(&Domain);

	Domain.SetHasHHR(HasHHR != 0);
	//Domain.CheckParams();

	// Call the function defined in DriverMexProb.h
	ParseSolverParamsAndOptimizing(prhs[4], &Prob, &SphereX, nullptr, plhs);

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
