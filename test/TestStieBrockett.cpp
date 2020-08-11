
#include "test/TestStieBrockett.h"

using namespace ROPTLIB;

/*User-specified linesearch algorithm*/
double LinesearchInput(integer iter, const Variable &x1, const Vector &exeta1, realdp initialstepsize, realdp initialslope, const Problem *prob, const Solvers *solver)
{
    return 1;
}

/*User-specified stopping criterion*/
bool MyStop(const Variable &x, const Vector &funSeries, integer lengthSeries, realdp finalval, realdp initval, const Problem *prob, const Solvers *solver)
{
    return (finalval / initval < 1e-6);
};

void testStieBrockett(void)
{
	unsigned tt = (unsigned)time(NULL);
	tt = 2; /*The following test is only for random seed zero*/
    std::cout << "seed SB:" << tt << std::endl;//---
//    tt = 1586166957; //
//    tt = 1;
	genrandseed(tt);

	// size of the Stiefel manifold
	integer n = 5, p = 2;
	// Generate the matrices in the Brockett problem.
    Vector B(n, n), D(p);
    B.RandGaussian();
    B = B + B.GetTranspose();
    realdp *Dptr = D.ObtainWriteEntireData();
	/*D is a diagonal matrix.*/
	for (integer i = 0; i < p; i++)
		Dptr[i] = static_cast<realdp> (i + 1);


//	Variable StieX(n, p);
//	StieX.RandInManifold();

	// Define the manifold
	Stiefel Domain(n, p);
	//Grassmann Domain(n, p);
//    Domain.ChooseParamsSet5();
    Variable StieX = Domain.RandominManifold();

							// Define the Brockett problem
	StieBrockett Prob(B, D);
	/*The domain of the problem is a Stiefel manifold*/
	Prob.SetDomain(&Domain);
	//Domain.ChooseParamsSet5();
	/*Output the parameters of the domain manifold*/
    Domain.ChooseParamsSet4();
//    Domain.SetHasHHR(true); /*set whether the manifold uses the idea in [HGA2015, Section 4.3] or not*/
	Domain.CheckParams();
    

//	SphereTx DomainPH(&Domain, &StieX);
//	SphereTxRQ ProbHess(&Domain, &StieX, &Prob);
//	ProbHess.SetDomain(&DomainPH);
//	Variable TV0 = DomainPH.RandominManifold();
////	ProbHess.CheckGradHessian(TV0);
//	RTRNewton RTRNewtonsolver(&ProbHess, &TV0);
//	RTRNewtonsolver.Verbose = FINALRESULT;
//	RTRNewtonsolver.Run();
    
//    Domain.CheckIntrExtr(StieX);
//	Domain.CheckRetraction(StieX);
	Domain.CheckDiffRetraction(StieX, true);
//	Domain.CheckLockingCondition(StieX);
//	Domain.CheckcoTangentVector(StieX);
//	Domain.CheckIsometryofVectorTransport(StieX);
//	Domain.CheckIsometryofInvVectorTransport(StieX);
//	Domain.CheckVecTranComposeInverseVecTran(StieX);
//	Domain.CheckTranHInvTran(StieX);
//	Domain.CheckHaddScaledRank1OPE(StieX);
//    return;
//    Prob.CheckGradHessian(StieX);
    
    LRBFGS *RSDsolver = new LRBFGS(&Prob, &StieX);
    RSDsolver->Verbose = ITERRESULT;//--- FINALRESULT;
    RSDsolver->LineSearch_LS = LSSM_INPUTFUN;
    RSDsolver->LinesearchInput = &LinesearchInput;
    RSDsolver->IsPureLSInput = false;
    RSDsolver->StopPtr = &MyStop;
    RSDsolver->Max_Iteration = 100;
    RSDsolver->OutputGap = 1;
    RSDsolver->CheckParams();
    RSDsolver->Run();
    
//    Prob.CheckGradHessian(RSDsolver->GetXopt());
    
//    Vector xopt = RSDsolver->GetXopt();
//    Prob.MinMaxEigValHess(StieX).Print("EigValHess:");
//    Prob.MinMaxEigValHess(xopt).Print("EigValHess:");
    delete RSDsolver;
}

///*We don't have to a line search algorithm defined in the solvers. The line seach algorithm can be defined
//here:*/
//realdp StieBrockettLinesearchInput(integer iter, Variable *x1, Vector *eta1, realdp initialstepsize, realdp initialslope, const Problem *prob, const Solvers *solver)
//{ /*For example, simply use one to be the stepsize*/
//
//	const StieBrockett *P = dynamic_cast<StieBrockett *> (const_cast<Problem*> (prob));
//	const realdp *etaTV = eta1->ObtainReadData();
//	const realdp *xM = x1->ObtainReadData();
//
//	integer n = P->n, p = P->p, length = n * p;
//	realdp denor, nume;
//	realdp *B = P->B, *D = P->D;
//	realdp *tmp = new realdp [n * p];
//	gemm_(GLOBAL::N, GLOBAL::N, &n, &p, &n, &GLOBAL::DONE, B, &n, const_cast<realdp *> (etaTV), &n, &GLOBAL::DZERO, tmp, &n);
//
//	denor = dot_(&length, tmp, &GLOBAL::IONE, const_cast<realdp *> (etaTV), &GLOBAL::IONE);
//	nume = dot_(&length, tmp, &GLOBAL::IONE, const_cast<realdp *> (xM), &GLOBAL::IONE);
//
//	delete[] tmp;
//	return (- nume / denor < 0) ? 1 : - nume / denor;
//}

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
	realdp *B, *D, *X;
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

	genrandseed(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;

	// Obtain an initial iterate by taking the Q factor of qr decomposition
	Variable StieX(n, p);
	realdp *StieXptr = StieX.ObtainWriteEntireData();
	for (integer i = 0; i < n * p; i++)
		StieXptr[i] = X[i];

	// Define the manifold
	Stiefel Domain(n, p);
	if (Paramset == 1)
		Domain.ChooseParamsSet1();
	else if (Paramset == 2)
		Domain.ChooseParamsSet2();
	else if (Paramset == 3)
		Domain.ChooseParamsSet3();
    else if (Paramset == 4)
        Domain.ChooseParamsSet4();
    else if (Paramset == 5)
        Domain.ChooseParamsSet5();

    Vector BB(n, n), DD(p);
    realdp *BBptr = BB.ObtainWriteEntireData();
    for(integer i = 0; i < n * n; i++)
        BBptr[i] = B[i];
    realdp *DDptr = DD.ObtainWriteEntireData();
    for(integer i = 0; i < p; i++)
        DDptr[i] = D[i];
    
	// Define the Brockett problem
	StieBrockett Prob(BB, DD);
	Prob.SetDomain(&Domain);

	Domain.SetHasHHR(HasHHR != 0);
	//Domain.CheckParams();

	// Call the function defined in DriverMexProb.h
	ParseSolverParamsAndOptimizing(prhs[5], &Prob, &StieX, plhs);

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
