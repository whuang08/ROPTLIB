
#include "test/TestEucFrechetMean.h"

using namespace ROPTLIB;

void testEucFrechetMean(void)
{
	// size of the samples and number of samples
	integer num = 6;
	integer dim = 4;

	// Generate the matrices in the Euclidean Frechetmean problem.
	double *Weights, *Data;
	double sumW;
	Weights = new double[num + num * dim];
	Data = Weights + num;
	for (integer i = 0; i < num; i++)
		Weights[i] = genrandreal() * 1e2 + 1e-4;
	Weights[0] = 1e-2;
	Weights[1] = 1e4;
	sumW = 0;
	for (integer i = 0; i < num; i++)
		sumW += Weights[i];
	for (integer i = 0; i < num; i++)
		Weights[i] /= sumW;
	for (integer i = 0; i < num * dim; i++)
		Data[i] = genrandnormal();

	testEucFrechetMean(Data, Weights, num, dim);

	delete[] Weights;
}

void testEucFrechetMean(double *Data, double *Weights, integer num, integer dim, double *X, double *Xopt)
{
	// Obtain an initial iterate
	EucVariable EucX(dim);
	if (X == nullptr)
	{
		EucX.RandInManifold();
	}
	else
	{
		double *EucXptr = EucX.ObtainWriteEntireData();
		for (integer i = 0; i < dim; i++)
			EucXptr[i] = X[i];
	}

	// Define the manifold
	Euclidean Domain(dim);

	// Define the problem
	EucFrechetMean Prob(Weights, Data, num, dim);
	Prob.SetDomain(&Domain);

	//// test RSD
	//for (integer i = 0; i < INPUTFUN; i++)
	//{
	//	RSD *RSDsolver = new RSD(&Prob, &EucX);
	//	RSDsolver->LineSearch_LS = static_cast<LSAlgo> (i);
	//	RSDsolver->Debug = FINALRESULT;
	//	//RSDsolver->CheckParams();
	//	RSDsolver->Run();
	//	delete RSDsolver;
	//}
	//// test RNewton
	//for (integer i = 0; i < INPUTFUN; i++)
	//{
	//	RNewton *RNewtonsolver = new RNewton(&Prob, &EucX);
	//	RNewtonsolver->LineSearch_LS = static_cast<LSAlgo> (i);
	//	RNewtonsolver->Debug = FINALRESULT;
	//	//RNewtonsolver->CheckParams();
	//	RNewtonsolver->Run();
	//	delete RNewtonsolver;
	//}

	//// test RCG
	//for (integer i = 0; i < RCGMETHODSLENGTH; i++)
	//{
	//	RCG *RCGsolver = new RCG(&Prob, &EucX);
	//	RCGsolver->RCGmethod = static_cast<RCGmethods> (i);
	//	RCGsolver->LineSearch_LS = STRONGWOLFE;
	//	RCGsolver->LS_beta = 0.1;
	//	RCGsolver->Debug = FINALRESULT;
	//	//RCGsolver->CheckParams();
	//	RCGsolver->Run();
	//	delete RCGsolver;
	//}

	//// test RBroydenFamily
	//for (integer i = 0; i < INPUTFUN; i++)
	//{
	//	RBroydenFamily *RBroydenFamilysolver = new RBroydenFamily(&Prob, &EucX);
	//	RBroydenFamilysolver->LineSearch_LS = static_cast<LSAlgo> (i);
	//	RBroydenFamilysolver->Debug = FINALRESULT;
	//	//RBroydenFamilysolver->CheckParams();
	//	RBroydenFamilysolver->Run();
	//	delete RBroydenFamilysolver;
	//}

	//// test RWRBFGS
	//for (integer i = 0; i < INPUTFUN; i++)
	//{
	//	RWRBFGS *RWRBFGSsolver = new RWRBFGS(&Prob, &EucX);
	//	RWRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
	//	RWRBFGSsolver->Debug = FINALRESULT;
	//	//RWRBFGSsolver->CheckParams();
	//	RWRBFGSsolver->Run();
	//	delete RWRBFGSsolver;
	//}

	//// test RBFGS
	//for (integer i = 0; i < INPUTFUN; i++)
	//{
	//	RBFGS *RBFGSsolver = new RBFGS(&Prob, &EucX);
	//	RBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
	//	RBFGSsolver->Debug = FINALRESULT;
	//	//RBFGSsolver->CheckParams();
	//	RBFGSsolver->Run();
	//	delete RBFGSsolver;
	//}

	// test LRBFGS
	for (integer i = 0; i < 1; i++) // INPUTFUN
	{
		LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, &EucX);
		LRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		LRBFGSsolver->Debug = FINALRESULT;
		//LRBFGSsolver->CheckParams();
		LRBFGSsolver->Run();
		if (LRBFGSsolver->Getnormgfgf0() < 1e-10)
			printf("SUCCESS!\n");
		else
			printf("FAIL!\n");
		delete LRBFGSsolver;
	}

	//RTRSD RTRSDsolver(&Prob, &EucX);
	//RTRSDsolver.Debug = FINALRESULT;
	////RTRSDsolver.CheckParams();
	//RTRSDsolver.Run();

	//RTRNewton RTRNewtonsolver(&Prob, &EucX);
	//RTRNewtonsolver.Debug = FINALRESULT;
	////RTRNewtonsolver.CheckParams();
	//RTRNewtonsolver.Run();

	//RTRSR1 RTRSR1solver(&Prob, &EucX);
	//RTRSR1solver.Debug = FINALRESULT;
	////RTRSR1solver.CheckParams();
	//RTRSR1solver.Run();

	//LRTRSR1 LRTRSR1solver(&Prob, &EucX);
	//LRTRSR1solver.Debug = FINALRESULT;
	////LRTRSR1solver.CheckParams();
	//LRTRSR1solver.Run();

	//// Check gradient and Hessian
	//Prob.CheckGradHessian(&EucX);
	//const Variable *xopt = RTRNewtonsolver.GetXopt();
	//Prob.CheckGradHessian(xopt);

	//if (Xopt != nullptr)
	//{
	//	const double *xoptptr = xopt->ObtainReadData();
	//	for (integer i = 0; i < dim; i++)
	//		Xopt[i] = xoptptr[i];
	//}
};

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs < 3)
    {
        mexErrMsgTxt("The number of arguments should be at least three.\n");
    }
    
	double *Data, *Weights, *X, *Xopt;
	Data = mxGetPr(prhs[0]);
	Weights = mxGetPr(prhs[1]);
	X = mxGetPr(prhs[2]);
	/* dimensions of input matrices */
	integer dim, num;
	dim = mxGetM(prhs[0]);
	num = mxGetN(prhs[0]);
    if(mxGetM(prhs[2]) != dim || mxGetN(prhs[2]) != 1)
    {
        mexErrMsgTxt("The size of initial X is not correct!\n");
    }
    if(mxGetM(prhs[1]) != num || mxGetN(prhs[1]) != 1)
    {
        mexErrMsgTxt("The size of the weights is not correct!\n");
    }
    
	printf("(dim, num):%d, %d\n", dim, num);

	/*create output matrix*/
	plhs[0] = mxCreateDoubleMatrix(dim, 1, mxREAL);
	Xopt = mxGetPr(plhs[0]);

	genrandseed(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;
	testEucFrechetMean(Data, Weights, num, dim, X, Xopt);
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
