
#include "test/TestCStieBrockett.h"

using namespace ROPTLIB;

void testCStieBrockett(void)
{
	unsigned tt = (unsigned)time(NULL);
	tt = 2; /*The following test is only for random seed zero*/
    std::cout << "seed SB:" << tt << std::endl;//---
//    tt = 1586166957; //
//    tt = 1;
	genrandseed(tt);

	// size of the Stiefel manifold
	integer n = 3, p = 2;
	// Generate the matrices in the Brockett problem.
    Vector B(n, n, "complex"), D(p, "complex");
    B.RandGaussian();
    B = B + B.GetTranspose();
    realdpcomplex *Dptr = (realdpcomplex *) D.ObtainWriteEntireData();
	/*D is a diagonal matrix.*/
	for (integer i = 0; i < p; i++)
    {
        Dptr[i].r = static_cast<realdp> (i + 1);
        Dptr[i].i = static_cast<realdp> (0);
    }


//	Variable StieX(n, p);
//	StieX.RandInManifold();

	// Define the manifold
	CStiefel Domain(n, p);
	//Grassmann Domain(n, p);
//    Domain.ChooseParamsSet5();
    Variable CStieX = Domain.RandominManifold();

							// Define the Brockett problem
	CStieBrockett Prob(B, D);
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
    
//    Domain.CheckIntrExtr(CStieX);
//	Domain.CheckRetraction(CStieX);
//	Domain.CheckDiffRetraction(CStieX, true);
//	Domain.CheckLockingCondition(CStieX);
//	Domain.CheckcoTangentVector(CStieX);
//	Domain.CheckIsometryofVectorTransport(CStieX);
//	Domain.CheckIsometryofInvVectorTransport(CStieX);
//	Domain.CheckVecTranComposeInverseVecTran(CStieX);
//	Domain.CheckTranHInvTran(CStieX);
//	Domain.CheckHaddScaledRank1OPE(CStieX);
//    return;
    
//    Prob.CheckGradHessian(CStieX);
    
    LRBFGS *RSDsolver = new LRBFGS(&Prob, &CStieX);
    RSDsolver->Verbose = FINALRESULT;//--- FINALRESULT;
//        RSDsolver->InitSteptype = LSSM_BBSTEP;
//    RSDsolver->LengthSY = 2;
    RSDsolver->Max_Iteration = 100;
    RSDsolver->OutputGap = 1;
    RSDsolver->CheckParams();
    RSDsolver->Run();
    B.Print("B:");//---
//    Prob.CheckGradHessian(RSDsolver->GetXopt());
    RSDsolver->GetXopt().Print("xopt:");//--
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
	B = (realdp *) mxGetComplexDoubles(prhs[0]);
	D = (realdp *) mxGetComplexDoubles(prhs[1]);
	X = (realdp *) mxGetComplexDoubles(prhs[2]);
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
	Variable StieX(n, p, "complex");
	realdp *StieXptr = StieX.ObtainWriteEntireData();
	for (integer i = 0; i < n * p * 2; i++)
		StieXptr[i] = X[i];

	// Define the manifold
	CStiefel Domain(n, p);
	if (Paramset == 1)
		Domain.ChooseParamsSet1();
	else if (Paramset == 2)
		Domain.ChooseParamsSet2();
    else if (Paramset == 3)
        Domain.ChooseParamsSet3();
    else if (Paramset == 4)
        Domain.ChooseParamsSet4();

    Vector BB(n, n, "complex"), DD(p, "complex");
    realdp *BBptr = BB.ObtainWriteEntireData();
    for(integer i = 0; i < n * n * 2; i++)
        BBptr[i] = B[i];
    realdp *DDptr = DD.ObtainWriteEntireData();
    for(integer i = 0; i < p * 2; i++)
        DDptr[i] = D[i];
    
	// Define the Brockett problem
	CStieBrockett Prob(BB, DD);
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
