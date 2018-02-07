
#include "test/TestStieBrockett.h"

using namespace ROPTLIB;

void testStieBrockettMore(void)
{
	// size of the Stiefel manifold
	integer n = 10, p = 4;

	// Generate the matrices in the Brockett problem.
	double *B = new double[n * n + p];
	double *D = B + n * n;
	/*B is an n by n matrix*/
	for (integer i = 0; i < n; i++)
	{
		for (integer j = i; j < n; j++)
		{
			B[i + j * n] = genrandnormal();
			B[j + i * n] = B[i + j * n];
		}
	}
	/*D is a diagonal matrix.*/
	for (integer i = 0; i < p; i++)
		D[i] = static_cast<double> (i + 1);


	StieVariable StieX(n, p);
	StieX.RandInManifold();

	// Define the manifold
	Stiefel Domain(n, p);
	//Grassmann Domain(n, p);
	Domain.SetHasHHR(true); /*set whether the manifold uses the idea in [HGA2015, Section 4.3] or not*/

							// Define the Brockett problem
	StieBrockett Prob(B, D, n, p);
	/*The domain of the problem is a Stiefel manifold*/
	Prob.SetDomain(&Domain);

	/*Output the parameters of the domain manifold*/
	//Domain.CheckParams();


	//SphereTx DomainPH(&Domain, &StieX);
	//SphereTxRQ ProbHess(&Domain, &StieX, &Prob);
	//ProbHess.SetDomain(&DomainPH);
	//Variable *TV0 = DomainPH.RandominManifold();
	//ProbHess.CheckGradHessian(TV0);
	//RTRNewton RTRNewtonsolver(&ProbHess, TV0);
	//RTRNewtonsolver.Debug = FINALRESULT;
	//RTRNewtonsolver.Run();
	//delete TV0;




	//Domain.CheckRetraction(&StieX);
	//Domain.CheckDiffRetraction(&StieX);
	//Domain.CheckLockingCondition(&StieX);
	//Domain.CheckcoTangentVector(&StieX);
	//Domain.CheckIsometryofVectorTransport(&StieX);
	//Domain.CheckIsometryofInvVectorTransport(&StieX);
	//Domain.CheckVecTranComposeInverseVecTran(&StieX);
	//Domain.CheckTranHInvTran(&StieX);
	//Domain.CheckHaddScaledRank1OPE(&StieX);

	//RBFGSLPSub *RBFGSLPSubsolver = new RBFGSLPSub(&Prob, &StieX);
	//RBFGSLPSubsolver->Debug = FINALRESULT;
	//RBFGSLPSubsolver->OutputGap = 10;
	//RBFGSLPSubsolver->lambdaLower = 1e-3;
	//RBFGSLPSubsolver->lambdaUpper = 1e3;
	//RBFGSLPSubsolver->CheckParams();
	//RBFGSLPSubsolver->Run();
	//delete RBFGSLPSubsolver;
	//return;

	//RBFGSLPSubsolver = new RBFGSLPSub(&Prob, &StieX);
	//RBFGSLPSubsolver->Debug = FINALRESULT;
	//RBFGSLPSubsolver->OutputGap = 10;
	//RBFGSLPSubsolver->lambdaLower = 1e-7;
	//RBFGSLPSubsolver->lambdaUpper = 1e7;
	//RBFGSLPSubsolver->CheckParams();
	//RBFGSLPSubsolver->Run();
	//delete RBFGSLPSubsolver;


	//RSD *RSDsolver = new RSD(&Prob, &StieX);
	//RSDsolver->Debug = FINALRESULT;
	//RSDsolver->InitSteptype = ONESTEP;
	//RSDsolver->Max_Iteration = 2000;
	//RSDsolver->CheckParams();
	//RSDsolver->Run();
	//delete RSDsolver;

	//RBFGS *RBFGSsolver = new RBFGS(&Prob, &StieX);
	//RBFGSsolver->Debug = FINALRESULT;
	//RBFGSsolver->CheckParams();
	//RBFGSsolver->Run();
	//delete RBFGSsolver;

	//// test LRBFGS
	//printf("\n********************************Check all line search algorithm in LRBFGS*************************************\n");
	////Domain.ChooseStieParamsSet4();
	//for (integer i = 0; i < 1; i++)//INPUTFUN
	//{
	//	LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, &StieX);
	//	LRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
	//	LRBFGSsolver->Debug = FINALRESULT; //ITERRESULT;// 
	//	//LRBFGSsolver->OutputGap = 1;
	//	LRBFGSsolver->Max_Iteration = 500;
	//	//LRBFGSsolver->Accuracy = 1e-6;
	//	LRBFGSsolver->Tolerance = 1e-6;
	//	//LRBFGSsolver->Finalstepsize = 1;
	//	//LRBFGSsolver->Num_pre_funs = 0;
	//	//LRBFGSsolver->BBratio = 1;
	//	//LRBFGSsolver->Num_pre_BB = 0;
	//	//LRBFGSsolver->InitSteptype = ONESTEP;
	//	LRBFGSsolver->LengthSY = 4;
	//	//LRBFGSsolver->CheckParams();
	//	LRBFGSsolver->Run();

	//	//// Compute the smallest eigenvalue of the Hessian at root.
	//	//Variable *root = StieX.ConstructEmpty();
	//	//LRBFGSsolver->GetXopt()->CopyTo(root);
	//	//SphereTx DomainPH(&Domain, root);
	//	//SphereTxRQ ProbHess(&Domain, root, &Prob, true);
	//	//ProbHess.SetDomain(&DomainPH);
	//	//Variable *TV0 = DomainPH.RandominManifold();
	//	//RTRNewton RTRNewtonsolver(&ProbHess, TV0);
	//	//RTRNewtonsolver.Debug = FINALRESULT;
	//	//RTRNewtonsolver.Run();
	//	//delete root;
	//	//delete TV0;

	//	//// Check the correctness of gradient and Hessian at the initial iterate
	//	//Prob.CheckGradHessian(&StieX);
	//	//const Variable *xopt = LRBFGSsolver->GetXopt();
	//	//// Check the correctness of gradient and Hessian at the final iterate of RTRNewton method
	//	//Prob.CheckGradHessian(xopt);
	//	delete LRBFGSsolver;
	//}

	//// test RTRSD
	//printf("\n********************************Check RTRSD*************************************\n");
	//RTRSD RTRSDsolver(&Prob, &StieX);
	//RTRSDsolver.Debug = FINALRESULT;
	//RTRSDsolver.Max_Iteration = 500;
	////RTRSDsolver.CheckParams();
	//RTRSDsolver.Run();

	//// test RTRNewton
	//printf("\n********************************Check RTRNewton*************************************\n");
	//RTRNewton RTRNewtonsolver(&Prob, &StieX);
	//RTRNewtonsolver.Debug = FINALRESULT;
	//RTRNewtonsolver.Max_Iteration = 11;
	////RTRNewtonsolver.CheckParams();
	//RTRNewtonsolver.Run();

	//// test RTRSR1
	//printf("\n********************************Check RTRSR1*************************************\n");
	//RTRSR1 RTRSR1solver(&Prob, &StieX);
	//RTRSR1solver.Debug = FINALRESULT;
	//RTRSR1solver.Max_Iteration = 65;
	////RTRSR1solver.CheckParams();
	//RTRSR1solver.Run();

	// test LRTRSR1
	printf("\n********************************Check LRTRSR1*************************************\n");
	LRTRSR1 LRTRSR1solver(&Prob, &StieX);
	LRTRSR1solver.OutputGap = 1;
	LRTRSR1solver.Max_Iteration = 500;
	LRTRSR1solver.Debug = ITERRESULT;
	//LRTRSR1solver.CheckParams();
	LRTRSR1solver.Tolerance = 1e-6;
	LRTRSR1solver.Run();

	//// Check the correctness of gradient and Hessian at the initial iterate
	//Prob.CheckGradHessian(&StieX);
	//const Variable *xopt = RTRNewtonsolver.GetXopt();
	//// Check the correctness of gradient and Hessian at the final iterate of RTRNewton method
	//Prob.CheckGradHessian(xopt);

	////Output the optimizer obtained by RTRNewton method
	//if (Xopt != nullptr)
	//{
	//	const double *xoptptr = xopt->ObtainReadData();
	//	for (integer i = 0; i < n * p; i++)
	//		Xopt[i] = xoptptr[i];
	//}

	delete[] B;
}

void testStieBrockett(void)
{
	// size of the Stiefel manifold
	integer n = 10, p = 4;

	// Generate the matrices in the Brockett problem.
	double *B = new double[n * n + p];
	double *D = B + n * n;
	/*B is an n by n matrix*/
	for (integer i = 0; i < n; i++)
	{
		for (integer j = i; j < n; j++)
		{
			B[i + j * n] = genrandnormal();
			B[j + i * n] = B[i + j * n];
		}
	}
	/*D is a diagonal matrix.*/
	for (integer i = 0; i < p; i++)
		D[i] = static_cast<double> (i + 1);

	testStieBrockett(B, D, n, p);

	delete[] B;
}

/*We don't have to a line search algorithm defined in the solvers. The line seach algorithm can be defined 
here:*/
double StieBrockettLinesearchInput(integer iter, Variable *x1, Vector *eta1, double initialstepsize, double initialslope, const Problem *prob, const Solvers *solver)
{ /*For example, simply use one to be the stepsize*/

	const StieBrockett *P = dynamic_cast<StieBrockett *> (const_cast<Problem*> (prob));
	const double *etaTV = eta1->ObtainReadData();
	const double *xM = x1->ObtainReadData();

	integer n = P->n, p = P->p, length = n * p;
	double denor, nume;
	double *B = P->B, *D = P->D;
	double *tmp = new double [n * p];
	dgemm_(GLOBAL::N, GLOBAL::N, &n, &p, &n, &GLOBAL::DONE, B, &n, const_cast<double *> (etaTV), &n, &GLOBAL::DZERO, tmp, &n);

	denor = ddot_(&length, tmp, &GLOBAL::IONE, const_cast<double *> (etaTV), &GLOBAL::IONE);
	nume = ddot_(&length, tmp, &GLOBAL::IONE, const_cast<double *> (xM), &GLOBAL::IONE);

	delete[] tmp;
	return (- nume / denor < 0) ? 1 : - nume / denor;
}

void testStieBrockett(double *B, double *D, integer n, integer p, double *X, double *Xopt)
{
	StieVariable StieX(n, p);

	if (X == nullptr)
	{/*If X is not defined before, then obtain an initial iterate by taking the Q factor of qr decomposition*/
		StieX.RandInManifold();
	}
	else
	{/*Otherwise, using the input orthonormal matrix as the initial iterate*/
		double *StieXptr = StieX.ObtainWriteEntireData();
		for (integer i = 0; i < n * p; i++)
			StieXptr[i] = X[i];
	}

	// Define the manifold
	Stiefel Domain(n, p);
	//Grassmann Domain(n, p);
	Domain.SetHasHHR(true); /*set whether the manifold uses the idea in [HGA2015, Section 4.3] or not*/

	// Define the Brockett problem
	StieBrockett Prob(B, D, n, p);
	/*The domain of the problem is a Stiefel manifold*/
	Prob.SetDomain(&Domain);

	// test RSD
	integer RSDMI[5] = { 2000, 2000, 2000, 2000, 2000 };
	printf("********************************Check all line search algorithm in RSD*****************************************\n");
	for (integer i = 0; i < INPUTFUN; i++)
	{
		RSD *RSDsolver = new RSD(&Prob, &StieX);
		RSDsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		RSDsolver->Debug = FINALRESULT;
		RSDsolver->Max_Iteration = RSDMI[i];
		//RSDsolver->CheckParams();
		RSDsolver->Run();
		if (RSDsolver->Getnormgfgf0() < 1e-6)
			printf("SUCCESS!\n");
		else
			printf("FAIL!\n");
		delete RSDsolver;
	}

	// test RNewton
	integer RNewtonMI[5] = { 50, 50, 50, 50, 50 };
	printf("\n********************************Check all line search algorithm in RNewton*************************************\n");
	for (integer i = 0; i < INPUTFUN; i++) //LSALGOLENGTH
	{
		RNewton *RNewtonsolver = new RNewton(&Prob, &StieX);
		RNewtonsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		RNewtonsolver->Debug = FINALRESULT;
		/*Uncomment following two lines to use the linesearch algorithm defined by the function "LinesearchInput".*/
		//RNewtonsolver->LineSearch_LS = INPUTFUN;
		//RNewtonsolver->LinesearchInput = &LinesearchInput;
		RNewtonsolver->Max_Iteration = RNewtonMI[i];
		//RNewtonsolver->CheckParams();
		RNewtonsolver->Run();

		if (RNewtonsolver->Getnormgfgf0() < 1e-6)
			printf("SUCCESS!\n");
		else
			printf("FAIL!\n");

		delete RNewtonsolver;
	}

	// test RCG
	integer RCGMI[6] = { 2000, 2000, 2000, 2000, 2000, 2000 };
	printf("\n********************************Check all Formulas in RCG*************************************\n");
	for (integer i = 0; i < RCGMETHODSLENGTH; i++)
	{
		RCG *RCGsolver = new RCG(&Prob, &StieX);
		RCGsolver->RCGmethod = static_cast<RCGmethods> (i);
		RCGsolver->LineSearch_LS = ARMIJO;
		//RCGsolver->InitSteptype = QUADINTMOD;
		RCGsolver->Debug = FINALRESULT;
		RCGsolver->Max_Iteration = RCGMI[i];
		//RCGsolver->CheckParams();
		RCGsolver->Run();
		if (RCGsolver->Getnormgfgf0() < 1e-6)
			printf("SUCCESS!\n");
		else
			printf("FAIL!\n");
		delete RCGsolver;
	}

	// test RBroydenFamily
	integer RBroydenFamilyMI[5] = { 200, 200, 200, 200, 200 };
	printf("\n********************************Check all Formulas in RBroydenFamily*************************************\n");
	for (integer i = 0; i < INPUTFUN; i++)
	{
		RBroydenFamily *RBroydenFamilysolver = new RBroydenFamily(&Prob, &StieX);
		RBroydenFamilysolver->LineSearch_LS = static_cast<LSAlgo> (i);
		RBroydenFamilysolver->Debug = FINALRESULT;
		RBroydenFamilysolver->Max_Iteration = RBroydenFamilyMI[i];
		//RBroydenFamilysolver->CheckParams();
		RBroydenFamilysolver->Run();
		if (RBroydenFamilysolver->Getnormgfgf0() < 1e-6)
			printf("SUCCESS!\n");
		else
			printf("FAIL!\n");
		delete RBroydenFamilysolver;
	}

	// test RWRBFGS
	integer RWRBFGSMI[5] = { 200, 200, 200, 200, 200 };
	printf("\n********************************Check all line search algorithm in RWRBFGS*************************************\n");
	for (integer i = 0; i < INPUTFUN; i++)
	{
		RWRBFGS *RWRBFGSsolver = new RWRBFGS(&Prob, &StieX);
		RWRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		RWRBFGSsolver->Debug = FINALRESULT; //ITERRESULT;//
		RWRBFGSsolver->Max_Iteration = RWRBFGSMI[i];
		//RWRBFGSsolver->CheckParams();
		RWRBFGSsolver->Run();
		if (RWRBFGSsolver->Getnormgfgf0() < 1e-6)
			printf("SUCCESS!\n");
		else
			printf("FAIL!\n");
		delete RWRBFGSsolver;
	}

	// test RBFGS
	integer RBFGSMI[5] = { 200, 200, 200, 200, 200 };
	printf("\n********************************Check all line search algorithm in RBFGS*************************************\n");
	for (integer i = 0; i < INPUTFUN; i++)
	{
		RBFGS *RBFGSsolver = new RBFGS(&Prob, &StieX);
		RBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		RBFGSsolver->Debug = FINALRESULT;
		RBFGSsolver->Max_Iteration = RBFGSMI[i];
		//RBFGSsolver->CheckParams();
		RBFGSsolver->Run();
		if (RBFGSsolver->Getnormgfgf0() < 1e-6)
			printf("SUCCESS!\n");
		else
			printf("FAIL!\n");
		delete RBFGSsolver;
	}

	// test LRBFGS
	integer LRBFGSMI[5] = { 200, 200, 200, 200, 200 };
	printf("\n********************************Check all line search algorithm in LRBFGS*************************************\n");
	//Domain.ChooseStieParamsSet4();
	for (integer i = 0; i < INPUTFUN; i++)//
	{
		LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, &StieX);
		LRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		LRBFGSsolver->Debug = FINALRESULT; //ITERRESULT;// 
		LRBFGSsolver->OutputGap = 1;
		LRBFGSsolver->Max_Iteration = LRBFGSMI[i];
		LRBFGSsolver->Accuracy = 1e-6;
		LRBFGSsolver->Tolerance = 1e-6;
		LRBFGSsolver->Finalstepsize = 1;
		LRBFGSsolver->Num_pre_funs = 0;
		LRBFGSsolver->BBratio = 1;
		LRBFGSsolver->Num_pre_BB = 0;
		LRBFGSsolver->InitSteptype = ONESTEP;
		LRBFGSsolver->LengthSY = 4;
		//LRBFGSsolver->CheckParams();
		LRBFGSsolver->Run();

		if (LRBFGSsolver->Getnormgfgf0() < 1e-6)
			printf("SUCCESS!\n");
		else
			printf("FAIL!\n");

		delete LRBFGSsolver;
	}

	// test RTRSD
	printf("\n********************************Check RTRSD*************************************\n");
	RTRSD RTRSDsolver(&Prob, &StieX);
	RTRSDsolver.Debug = FINALRESULT;
	RTRSDsolver.Max_Iteration = 2000;
	//RTRSDsolver.CheckParams();
	RTRSDsolver.Run();
	if (RTRSDsolver.Getnormgfgf0() < 1e-6)
		printf("SUCCESS!\n");
	else
		printf("FAIL!\n");

	// test RTRNewton
	printf("\n********************************Check RTRNewton*************************************\n");
	RTRNewton RTRNewtonsolver(&Prob, &StieX);
	RTRNewtonsolver.Debug = FINALRESULT;
	RTRNewtonsolver.Max_Iteration = 20;
	//RTRNewtonsolver.CheckParams();
	RTRNewtonsolver.Run();
	if (RTRNewtonsolver.Getnormgfgf0() < 1e-6)
		printf("SUCCESS!\n");
	else
		printf("FAIL!\n");

	// test RTRSR1
	printf("\n********************************Check RTRSR1*************************************\n");
	RTRSR1 RTRSR1solver(&Prob, &StieX);
	RTRSR1solver.Debug = FINALRESULT;
	RTRSR1solver.Max_Iteration = 200;
	//RTRSR1solver.CheckParams();
	RTRSR1solver.Run();
	if (RTRSR1solver.Getnormgfgf0() < 1e-6)
		printf("SUCCESS!\n");
	else
		printf("FAIL!\n");

	// test LRTRSR1
	printf("\n********************************Check LRTRSR1*************************************\n");
	LRTRSR1 LRTRSR1solver(&Prob, &StieX);
	LRTRSR1solver.OutputGap = 1;
	LRTRSR1solver.Max_Iteration = 500;
	//LRTRSR1solver.Shrinked_tau = 0.1;
	//LRTRSR1solver.LengthSY = 4;
	LRTRSR1solver.Debug = FINALRESULT;//--- FINALRESULT;
	//LRTRSR1solver.CheckParams();
	LRTRSR1solver.Tolerance = 1e-6;
	LRTRSR1solver.Run();
	if (LRTRSR1solver.Getnormgfgf0() < 1e-6)
		printf("SUCCESS!\n");
	else
		printf("FAIL!\n");
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
	if(nrhs < 6)
	{
		mexErrMsgTxt("The number of arguments should be at least six.\n");
	}
	double *B, *D, *X, *Xopt;
	B = mxGetPr(prhs[0]);
	D = mxGetPr(prhs[1]);
	X = mxGetPr(prhs[2]);
	/* dimensions of input matrices */
	integer p, n, HasHHR, Paramset;
	n = mxGetM(prhs[0]);
	p = mxGetM(prhs[1]);

	/*Check the correctness of the inputs*/
	if(mxGetN(prhs[0]) != n)
	{
		mexErrMsgTxt("The size of matrix B is not correct.\n");
	}
	if(mxGetN(prhs[1]) != 1)
	{
		mexErrMsgTxt("The size of the D is not correct!\n");
	}
	if(mxGetM(prhs[2]) != n || mxGetN(prhs[2]) != p)
	{
		mexErrMsgTxt("The size of the initial X is not correct!\n");
	}
	HasHHR = static_cast<integer> (mxGetScalar(prhs[3]));
	Paramset = static_cast<integer> (mxGetScalar(prhs[4]));

	printf("(n, p):%d,%d\n", n, p);

	/*create output matrix*/
	plhs[0] = mxCreateDoubleMatrix(n, p, mxREAL);
	Xopt = mxGetPr(plhs[0]);

	genrandseed(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;
	//	testStieBrockett(B, D, n, p, X, Xopt);

	// Obtain an initial iterate by taking the Q factor of qr decomposition
	StieVariable StieX(n, p);
	double *StieXptr = StieX.ObtainWriteEntireData();
	for (integer i = 0; i < n * p; i++)
		StieXptr[i] = X[i];

	StieVariable *Stiesoln = nullptr;
	if (nrhs >= 7)
	{
		double *soln = mxGetPr(prhs[6]);
		Stiesoln = new StieVariable(n, p);
		double *Stiesolnptr = Stiesoln->ObtainWriteEntireData();
		for (integer i = 0; i < n * p; i++)
		{
			Stiesolnptr[i] = soln[i];
		}
	}
	// Define the manifold
	Stiefel Domain(n, p);
	if (Paramset == 1)
		Domain.ChooseStieParamsSet1();
	else if (Paramset == 2)
		Domain.ChooseStieParamsSet2();
	else if (Paramset == 3)
		Domain.ChooseStieParamsSet3();

	// Define the Brockett problem
	StieBrockett Prob(B, D, n, p);
	Prob.SetDomain(&Domain);

	Domain.SetHasHHR(HasHHR != 0);
	//Domain.CheckParams();

	// Call the function defined in DriverMexProb.h
	ParseSolverParamsAndOptimizing(prhs[5], &Prob, &StieX, Stiesoln, plhs);

	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	delete Stiesoln;
	return;
}

#endif
