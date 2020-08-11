
#include "test/TestStieSPCA.h"

using namespace ROPTLIB;

void testStieSPCA(void)
{
	// size of the Stiefel manifold
    integer n = 5, m = 4, p = 3;
//    integer n = 2000, m = 50, p = 5;
    realdp lambda = 1;

	// Generate the matrices in the Brockett problem.
    Vector B(m, n);
    realdp *Bptr = B.ObtainWriteEntireData();
	/*B is an n by n matrix*/
	for (integer i = 0; i < n * m; i++)
	{
        Bptr[i] = genrandnormal();
	}
    for(integer i = 0; i < n; i++)
    {
        realdp s = 0;
        for(integer j = 0; j < m; j++)
        {
            s += Bptr[j + i * m];
        }
        s /= m;
        for(integer j = 0; j < m; j++)
            Bptr[j + i * m] -= s;
        s = 0;
        for(integer j = 0; j < m; j++)
        {
            s += Bptr[j + i * m] * Bptr[j + i * m];
        }
        s = std::sqrt(s);
        for(integer j = 0; j < m; j++)
            Bptr[j + i * m] /= s;
    }
    
//    Variable StieX(n, p);
//
//    StieX.RandInManifold();
//    ForDebug::Print("B", B, m, n);//---
//    StieX.Print("StieX:");//---
    // Define the manifold
    Stiefel Domain(n, p);
    Domain.ChooseParamsSet5();
    Variable StieX = Domain.RandominManifold();

    integer lengthW = 1;
    // Define the SPCA problem
    StieSPCA Prob(B, lambda, n, m, p, lengthW); // DIAGRHESS // LIPSCHITZ
    /*The domain of the problem is a Stiefel manifold*/
    Prob.SetDomain(&Domain);
    Domain.CheckParams();
    
//    ARPG *ARPGsolver = new ARPG(&Prob, &StieX);
//    ARPGsolver->Max_Iteration = 5000;
//    ARPGsolver->Variant = REGULAR; //-- REGULAR; //--ADALIPSCHITZ;
//    ARPGsolver->CheckParams();
//    ARPGsolver->Run();
//    delete ARPGsolver;
    
    AManPG *AManPGsolver = new AManPG(&Prob, &StieX);
    AManPGsolver->Max_Iteration = 100;
    AManPGsolver->Verbose = ITERRESULT;
    AManPGsolver->Variant = LSPG_REGULAR; //-- LSPG_REGULAR; //-- LSPG_ADALIPSCHITZ;
    AManPGsolver->CheckParams();
    AManPGsolver->Run();
    delete AManPGsolver;
    
    ManPG *ManPGsolver = new ManPG(&Prob, &StieX);
    ManPGsolver->Max_Iteration = 100;
    ManPGsolver->Verbose = ITERRESULT;
    ManPGsolver->Variant = LSPG_ADALIPSCHITZ; //-- LSPG_REGULAR; //-- LSPG_ADALIPSCHITZ;
    ManPGsolver->CheckParams();
    ManPGsolver->Run();
    delete ManPGsolver;
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
	realdp *A, *X;
    realdp lambda;
	A = mxGetPr(prhs[0]);
    lambda = static_cast<realdp> (mxGetScalar(prhs[1]));
	X = mxGetPr(prhs[2]);
	/* dimensions of input matrices */
	integer m, n, p, HasHHR;
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);
    p = mxGetN(prhs[2]);

	/*Check the correctness of the inputs*/
	if(mxGetM(prhs[2]) != n)
	{
		mexErrMsgTxt("The size of matrix Xinit is not correct.\n");
	}
	HasHHR = 0;

	genrandseed(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;

    Vector AA(m, n);
    realdp *AAptr = AA.ObtainWriteEntireData();
    for(integer i = 0; i < m * n; i++)
        AAptr[i] = A[i];
    
    // Define the manifold
    Stiefel Domain(n, p);
    Domain.ChooseParamsSet5();

	// Obtain an initial iterate by taking the Q factor of qr decomposition
	Variable initX = Domain.RandominManifold();
	realdp *initXptr = initX.ObtainWriteEntireData();
	for (integer i = 0; i < n * p; i++)
		initXptr[i] = X[i];
    
    integer lengthW = 1;
	// Define the Brockett problem
    StieSPCA Prob(AA, lambda, n, m, p, lengthW); // DIAGRHESS // LIPSCHITZ
	Prob.SetDomain(&Domain);

	Domain.SetHasHHR(HasHHR != 0);
	//Domain.CheckParams();

	// Call the function defined in DriverMexProb.h
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
