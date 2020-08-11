#include "test/TestFRankQ2FMatCompletion.h"

using namespace ROPTLIB;

//void testLRMatrixCompletionMore(void)
//{
//	integer m = 10, n = 10, r = 2;
//
//	LowRank Domain(m, n, r);
//	//Domain.SetHasHHR(true);
//	LowRankVariable InitialX(m, n, r);
//	InitialX.RandInManifold();
//	//Domain.CheckParams();
//	//Domain.CheckIntrExtr(&InitialX);
//	//Domain.CheckRetraction(&InitialX);
//	//Domain.CheckcoTangentVector(&InitialX);
//	//Domain.CheckDiffRetraction(&InitialX, false);
//	//Domain.CheckIsometryofVectorTransport(&InitialX);
//
//	//Domain.CheckLockingCondition(&InitialX);
//	//Domain.CheckIsometryofInvVectorTransport(&InitialX);
//	//Domain.CheckVecTranComposeInverseVecTran(&InitialX);
//	//Domain.CheckTranHInvTran(&InitialX);
//	//return;
//
//	//InitialX.Print("initialX:");
//
//	// Generate the matrices in the matrix completion approximation problem.
//	integer dim = (m + n - r) * r;
//	integer nz = 2 * dim;
//	integer *ir = new integer[nz * 2];
//	integer *jc = ir + nz;
//
//	integer *tmpforidx = new integer[m * n];
//	for (integer i = 0; i < m * n; i++)
//		tmpforidx[i] = i;
//	/*nz number of indices*/
//	integer idx = 0, itmp;
//	for (integer i = 0; i < nz; i++)
//	{
//		/*idx is an integer in [0, m - i - 1]*/
//		idx = static_cast<integer> ((m * n - i) * genrandreal());
//		while (idx >= m * n - i)
//			idx = static_cast<integer> ((m * n - i) * genrandreal());
//		/*the chosen idx is put at the end of the array*/
//		itmp = tmpforidx[m * n - i - 1];
//		tmpforidx[m * n - i - 1] = tmpforidx[idx];
//		tmpforidx[idx] = itmp;
//	}
//	for (integer i = 0; i < nz; i++)
//	{
//		/*tmpforidx[nz - 1 - i]*/
//		ir[i] = static_cast<integer> (tmpforidx[nz - 1 - i] / n);
//		jc[i] = tmpforidx[nz - 1 - i] - n * ir[i];
//	}
//	delete[] tmpforidx;
//
//	integer mn = m * n, mr = m * r, nr = n * r;
//	realdp *A_U = new realdp[mr];
//	realdp *A_V = new realdp[nr];
//	for (integer i = 0; i < m * r; i++)
//	{
//		A_U[i] = genrandnormal();
//	}
//	for (integer i = 0; i < n * r; i++)
//	{
//		A_V[i] = genrandnormal();
//	}
//	realdp *V = new realdp[nz];
//	for (integer i = 0; i < nz; i++)
//	{
//		V[i] = 0;
//		for (integer j = 0; j < r; j++)
//		{
//			V[i] += A_U[ir[i] + j * m] * A_V[jc[i] + j * n];
//		}
//	}
//	delete[]A_U;
//	delete[]A_V;
//
//	LRMatrixCompletion Prob(ir, jc, V, nz, m, n, r);
//	Prob.SetDomain(&Domain);
//
//	//Prob.f(&InitialX);
//	//Vector *gf = Domain.GetEMPTYINTR()->ConstructEmpty();//---
//	//Prob.Grad(&InitialX, gf);//---
//	//delete gf;
//	//	Prob.CheckGradHessian(&InitialX);
//
//	//	//RSD *RSDsolver = new RSD(&Prob, &InitialX);
//	//	//RTRNewton *RSDsolver = new RTRNewton(&Prob, &InitialX);
//	//	RCG *RSDsolver = new RCG(&Prob, &InitialX);
//	//	//->LineSearch_LS = ARMIJO;
//	//	//RSDsolver->LS_beta = 0.01;
//	//	//RSDsolver->RCGmethod = DAI_YUAN;
//	//	RSDsolver->Debug = ITERRESULT;
//	//	RSDsolver->OutputGap = 100;
//	//	RSDsolver->Max_Iteration = 500;
//	////	RSDsolver->LineSearch_LS = EXACT;
//	//	RSDsolver->CheckParams();
//	//	//RSDsolver->Accuracy = 1e-6;
//	//	RSDsolver->Tolerance = 1e-10;
//	//	/*Uncomment following two lines to use the linesearch algorithm defined by the function "LinesearchInput".*/
//	//	RSDsolver->LineSearch_LS = INPUTFUN;
//	//	RSDsolver->LinesearchInput = &LinesearchInput;
//	//	RSDsolver->Run();
//	//	//Prob.CheckGradHessian(&InitialX);//--
//	////	Prob.CheckGradHessian(RSDsolver->GetXopt());//--
//	//
//	//	delete RSDsolver;
//
//	// test LRBFGS
//	//printf("********************************Check LRBFGS*************************************\n");
//	LRBFGS LRBFGSsolver(&Prob, &InitialX);
//	LRBFGSsolver.OutputGap = 1;
//	LRBFGSsolver.Max_Iteration = 500;
//	LRBFGSsolver.Debug = FINALRESULT;//--- FINALRESULT;
//									 //LRBFGSsolver.CheckParams();
//	LRBFGSsolver.Tolerance = static_cast<realdp> (1e-6);
//	LRBFGSsolver.Run();
//
//	// test LRTRSR1
//	printf("********************************Check LRTRSR1*************************************\n");
//	LRTRSR1 LRTRSR1solver(&Prob, &InitialX);
//	LRTRSR1solver.OutputGap = 1;
//	LRTRSR1solver.Max_Iteration = 500;
//	//LRTRSR1solver.Shrinked_tau = 0.1;
//	//LRTRSR1solver.LengthSY = 8;
//	LRTRSR1solver.Debug = FINALRESULT;//--- FINALRESULT;
//	LRTRSR1solver.Tolerance = static_cast<realdp> (1e-6);
//	//LRTRSR1solver.CheckParams();
//	LRTRSR1solver.Run();
//
//	delete[] V;
//	delete[] ir;
//};

void testFRankQ2FMatCompletion(void)
{
	integer m = 8, n = 7, r = 2;

	FixedRankQ2F Domain(m, n, r);
	//Domain.SetHasHHR(true);
	Variable InitialX = Domain.RandominManifold();

	// Generate the matrices in the matrix completion approximation problem.
	integer dim = (m + n - r) * r;
	integer nz = 2 * dim;
	integer *ir = new integer[nz * 2];
	integer *jc = ir + nz;

	integer *tmpforidx = new integer[m * n];
	for (integer i = 0; i < m * n; i++)
		tmpforidx[i] = i;
	/*nz number of indices*/
	integer idx = 0, itmp;
	for (integer i = 0; i < nz; i++)
	{
		/*idx is an integer in [0, m - i - 1]*/
		idx = static_cast<integer> ((m * n - i) * genrandreal());
		while (idx >= m * n - i)
			idx = static_cast<integer> ((m * n - i) * genrandreal());
		/*the chosen idx is put at the end of the array*/
		itmp = tmpforidx[m * n - i - 1];
		tmpforidx[m * n - i - 1] = tmpforidx[idx];
		tmpforidx[idx] = itmp;
	}
	for (integer i = 0; i < nz; i++)
	{
		/*tmpforidx[nz - 1 - i]*/
		ir[i] = static_cast<integer> (tmpforidx[nz - 1 - i] / n);
		jc[i] = tmpforidx[nz - 1 - i] - n * ir[i];
	}
	delete[] tmpforidx;

	integer mr = m * r, nr = n * r;
	realdp *A_U = new realdp[mr];
	realdp *A_V = new realdp[nr];
	for (integer i = 0; i < m * r; i++)
	{
		A_U[i] = genrandnormal();
	}
	for (integer i = 0; i < n * r; i++)
	{
		A_V[i] = genrandnormal();
	}
	realdp *V = new realdp[nz];
	for (integer i = 0; i < nz; i++)
	{
		V[i] = 0;
		for (integer j = 0; j < r; j++)
		{
			V[i] += A_U[ir[i] + j * m] * A_V[jc[i] + j * n];
		}
	}
	delete[]A_U;
	delete[]A_V;

    FRankQ2FMatCompletion Prob(ir, jc, V, nz, m, n, r);
	Prob.SetDomain(&Domain);
//    Prob.SetNumGradHess(true);

    Prob.CheckGradHessian(InitialX);
//    return;//---
	// test LRBFGS
	//printf("********************************Check LRBFGS*************************************\n");
    RTRNewton LRBFGSsolver(&Prob, &InitialX);
//    LRBFGS LRBFGSsolver(&Prob, &InitialX);
	LRBFGSsolver.OutputGap = 1;
	LRBFGSsolver.Max_Iteration = 500;
	LRBFGSsolver.Verbose = FINALRESULT;//--- FINALRESULT;
									 //LRBFGSsolver.CheckParams();
	LRBFGSsolver.Tolerance = static_cast<realdp> (1e-6);
	LRBFGSsolver.Run();

    Prob.MinMaxEigValHess(LRBFGSsolver.GetXopt()).Print("eigs:");//---
//    Prob.CheckGradHessian(LRBFGSsolver.GetXopt());

	if (LRBFGSsolver.Getnormgfgf0() < 1e-6)
		printf("SUCCESS!\n");
	else
		printf("FAIL!\n");

	delete[] V;
	delete[] ir;
};

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 5)
	{
		mexErrMsgTxt("The number of arguments should be at least five.\n");
	}
	realdp *A, *X, *Xopt;
	A = mxGetPr(prhs[0]);
	X = mxGetPr(prhs[1]);
	/* dimensions of input matrices */
	integer m, n, HasHHR, nzmax, r;
	size_t *ir, *jc;
	nzmax = mxGetNzmax(prhs[0]);
	ir = mxGetIr(prhs[0]);
	jc = mxGetJc(prhs[0]);
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);
	r = static_cast<integer> (mxGetScalar(prhs[2]));

	/*Check the correctness of the inputs*/
	if (mxGetM(prhs[1]) != (m + n) * r || mxGetN(prhs[1]) != 1)
	{
		mexErrMsgTxt("The size of the initial X is not correct!\n");
	}
	HasHHR = static_cast<integer> (mxGetScalar(prhs[3]));

	genrandseed(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;

	// Define the manifold
	FixedRankQ2F Domain(m, n, r);
    
    Variable initX = Domain.RandominManifold();
    realdp *initXptr = initX.ObtainWriteEntireData();
    for(integer i = 0; i < initX.Getlength(); i++)
    {
        initXptr[i] = X[i];
    }
    
	// Define the matrix completion problem
	integer *inir = new integer[2 * nzmax];
	integer *injc = inir + nzmax;

	for (integer i = 0; i < n; i++)
	{
		for (unsigned long long j = jc[i]; j < jc[i + 1]; j++)
		{
			/*row: ir[j], column: i, entry: A[j]*/
			inir[j] = ir[j];
			injc[j] = i;
		}
	}
    
    FRankQ2FMatCompletion Prob(inir, injc, A, nzmax, m, n, r);
    Prob.SetDomain(&Domain);

	Domain.SetHasHHR(HasHHR != 0);

	// Call the function defined in DriverMexProb.h
	ParseSolverParamsAndOptimizing(prhs[4], &Prob, &initX, plhs);

	delete[] inir;

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
