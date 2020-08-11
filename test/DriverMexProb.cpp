
#include "test/DriverMexProb.h"

#ifdef MATLAB_MEX_FILE

using namespace ROPTLIB;

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	CheckMemoryDeleted = new std::map<integer *, integer>;
//    std::cout << "tt1:" << std::endl;//---
//    printf("tt2:\n");//---
	DriverMexProb(nlhs, plhs, nrhs, prhs);
//    printf("h17\n");//--
	// check memory
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
//    printf("h18\n");//--
	return;
}

void DriverMexProb(int &nlhs, mxArray ** &plhs, int &nrhs, const mxArray ** &prhs)
{
//    printf("h1\n");//--
	if (nrhs < 6)
	{
		mexErrMsgTxt("The number of arguments should be at least six.\n");
	}
//    printf("h2\n");//--
	// Argument Checking:
	// First to third arguments should be function handles
	if (!mxIsClass(prhs[0], "function_handle"))
	{
		mexErrMsgTxt("The first input argument is not a function handle.");
	}
	// fifth to sixth arguments are structures
	// The fifth one is SolverParams
	// The sixth one is ManiParams
	if (!mxIsStruct(prhs[4]) || !mxIsStruct(prhs[5]))
	{
		mexErrMsgTxt("At least one of fifth to sixth input arguments is not a structure.");
	}
    
//    printf("h3\n");//--
	// Obtain manifold and iterate structure
	Manifold *domain = nullptr, **manifolds = nullptr;
	Variable initialX;
//	Element **elements;
	integer *powsinterval, numoftype, numoftotal;
    
//    printf("h4\n");//--
	if (!ParseManiParams(prhs[5], manifolds, numoftype, powsinterval))
	{
		mexErrMsgTxt("Parsing ManiParams fails.");
	}
//    printf("h5\n");//--
    numoftotal = powsinterval[numoftype];
    if(numoftotal > 1)
        domain = new MultiManifolds(manifolds, numoftype, powsinterval);
    else
        domain = manifolds[0];
    
//    printf("h6\n");//--
	mxArray *tmp = mexProblem::GetFieldbyName(prhs[5], 0, "IsCheckParams");
	if (tmp != nullptr)
	{
		if (fabs(mxGetScalar(tmp)) > std::numeric_limits<realdp>::epsilon()) // if the value is nonzero
		{
			domain->CheckParams();
		}
	}
//    printf("h7\n");//--
	bool HasHHR = false;
	if (nrhs >= 7)
	{
		if (!mxIsDouble(prhs[6]))
		{
			mexErrMsgTxt("Seventh input argument is not a scalar.");
		}

		HasHHR = (static_cast<integer> (mxGetScalar(prhs[6])) != 0);
	}
	domain->SetHasHHR(HasHHR);
    initialX = domain->RandominManifold();
//    printf("h8\n");//--
    
	// initialize the initial iterate
	if (nrhs >= 8)
	{
		if (!mxIsStruct(prhs[7]))
		{
			mexErrMsgTxt("Eighth input argument is not a structure.");
		}
		mexProblem::ObtainElementFromMxArray(&initialX, prhs[7]);
	}
//    printf("h9\n");//--

//    initialX.Print("initialX:");
	// 	initialX->Print("initialX", false);

	// Define the problem
	Problem *Prob = new mexProblem(prhs[0], prhs[1], prhs[2], prhs[3]);
//    printf("h10\n");//--
	Prob->SetDomain(domain);
//    Prob->CheckGradHessian(initialX);
	//	Vector *egf = initialX->ConstructEmpty();
	//	Prob->EucGrad(initialX, egf);
	//	delete egf;
	// solve the optimization problem
	ParseSolverParamsAndOptimizing(prhs[4], Prob, &initialX, plhs);
    
//    printf("h11\n");//--
	delete Prob;
//    printf("h12\n");//--
	if (numoftotal > 1)
		delete domain;
    
//    printf("h13\n");//--
	for (integer i = 0; i < numoftype; i++)
	{
		delete manifolds[i];
	}
//    printf("h14\n");//--
	delete[] manifolds;
//    printf("h15\n");//--
	delete[] powsinterval;
//    printf("h16\n");//--
};

#endif // end of MATLAB_MEX_FILE
