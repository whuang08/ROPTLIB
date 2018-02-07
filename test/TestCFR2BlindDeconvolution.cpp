#include "test/TestCFR2BlindDeconvolution.h"

#ifdef ROPTLIB_WITH_FFTW

using namespace ROPTLIB;

//void testfft()
//{
//	int i, j, bw, bw2_1, size, size2_1, nrow, ncol;
//	int data_is_real;
//	int cutoff;
//	int rank, howmany_rank;
//	double *rresult, *iresult, *rdata, *idata;
//	double *workspace, *weights;
//
//	fftw_plan dctPlan;
//	fftw_plan fftPlan;
//	fftw_iodim dims[1], howmany_dims[1];
//
//	bw = 2;
//	weights = new double[4 * bw];// (double *)malloc(sizeof(double) * 4 * bw);
//	rdata = new double[5 * bw];// (double *)malloc(sizeof(double) * 5 * bw);
//	dctPlan = fftw_plan_r2r_1d(2 * bw, weights, rdata, FFTW_REDFT10, FFTW_ESTIMATE);
//	delete weights;
//	delete rdata;
//}

void testCFR2BlindDeconvolution(void)
{
	integer L = 16, K = 2, N = 4, r = 1;

	CFixedRank2Factors Domain(K, N, r);
	Domain.SetHasHHR(false);
	CFR2Variable InitialX(K, N, r);
	InitialX.RandInManifold();
	//Domain.CheckParams();

	//Domain.CheckIntrExtr(&InitialX);
	//Domain.CheckRetraction(&InitialX);
	//Domain.CheckcoTangentVector(&InitialX);
	//Domain.CheckDiffRetraction(&InitialX, false);
	//Domain.CheckIsometryofVectorTransport(&InitialX);

	//Domain.CheckLockingCondition(&InitialX);
	//Domain.CheckIsometryofInvVectorTransport(&InitialX);
	//Domain.CheckVecTranComposeInverseVecTran(&InitialX);
	//Domain.CheckTranHInvTran(&InitialX);




	// Generate the matrices in the Low rank approximation problem.
	double *y = new double[L * 2 + K * L * 2 + L * N * 2];
	double *B = y + L * 2;
	double *C = nullptr;//-- B + K * L * 2;
	for (integer i = 0; i < L * 2 + K * L * 2 + L * N * 2; i++)
		y[i] = genrandnormal();

	//ForDebug::Print("y:", y, L * 2);
	//ForDebug::Print("B:", B, L * 2, K);
	//ForDebug::Print("C:", C, L * 2, N);
	//InitialX.Print("InitX:");//---

	CFR2BlindDeconvolution Prob(y, B, 0, nullptr, nullptr, false, C, 0, nullptr, nullptr, 0, L, K, N, r, 0, 1, 1);
	Prob.SetDomain(&Domain);

	//Prob.CheckGradHessian(&InitialX);
	LRBFGS *RSDsolver = new LRBFGS(&Prob, &InitialX);
	RSDsolver->Debug = FINALRESULT;
	//RSDsolver->OutputGap = 100;
	RSDsolver->Max_Iteration = 50;
	//RSDsolver->CheckParams();
	RSDsolver->Accuracy = 1e-6;
	RSDsolver->Finalstepsize = 1;
	RSDsolver->Tolerance = 1e-6;
	RSDsolver->Run();
	//Prob.CheckGradHessian(&InitialX);
	//Prob.CheckGradHessian(RSDsolver->GetXopt());
	if (RSDsolver->Getnormgfgf0() < 1e-6)
		printf("SUCCESS!\n");
	else
		printf("FAIL!\n");
	delete RSDsolver;

	delete[] y;
};

void testCFR2BlindDeconvolutionSparse(void)
{
	//integer n = 8;
	//double *v = new double[n * 2];
	//for (integer i = 0; i < n * 2; i++)
	//	v[i] = i;
	//haarFWT_1d(n, (doublecomplex *)v); /*wavedec in Matlab*/

	//for (integer i = 0; i < n * 2; i++)
	//	std::cout << v[i] << std::endl;//---

	//haarFWT_1d_inverse(n, (doublecomplex *)v); /*waverec in Matlab*/

	//for (integer i = 0; i < n * 2; i++)
	//	std::cout << v[i] << std::endl;//---
	//delete[] v;
	//return;

	//integer n1 = 4, n2 = 4;
	//doublecomplex *vv = new doublecomplex[4 * 4];
	//for (integer i = 0; i < n1 * n2; i++)
	//{
	//	vv[i].r = i;
	//	vv[i].i = i;
	//}
	//std::cout << "original:" << std::endl;//---
	//for (integer i = 0; i < n1; i++)
	//{
	//	for (integer j = 0; j < n2; j++)
	//	{
	//		std::cout << vv[i + j * n1].r << "+ i " << vv[i + j * n1].i << "\t";
	//	}
	//	std::cout << std::endl;//---
	//}
	//haarFWT_2d(4, 4, vv);
	//std::cout << "Haar transform:" << std::endl;//---
	//for (integer i = 0; i < n1; i++)
	//{
	//	for (integer j = 0; j < n2; j++)
	//	{
	//		std::cout << vv[i + j * n1].r << "+ i " << vv[i + j * n1].i << "\t";
	//	}
	//	std::cout << std::endl;//---
	//}
	//std::cout << "Haar inverse transform:" << std::endl;//---
	//haarFWT_2d_inverse(n1, n2, vv);

	//for (integer i = 0; i < n1; i++)
	//{
	//	for (integer j = 0; j < n2; j++)
	//	{
	//		std::cout << vv[i + j * n1].r << "+ i " << vv[i + j * n1].i << "\t";
	//	}
	//	std::cout << std::endl;//---
	//}
	//delete[] vv;
	//return;

	integer K = 2, N = 4, r = 1, L = 16;//-- 3 * (K + N);

	CFixedRank2Factors Domain(K, N, r);
	Domain.SetHasHHR(false);
	CFR2Variable InitialX(K, N, r);
	InitialX.RandInManifold();
	//Domain.CheckParams();

	//Domain.CheckIntrExtr(&InitialX);
	//Domain.CheckRetraction(&InitialX);
	//Domain.CheckcoTangentVector(&InitialX);
	//Domain.CheckDiffRetraction(&InitialX, false);
	//Domain.CheckIsometryofVectorTransport(&InitialX);

	//Domain.CheckLockingCondition(&InitialX);
	//Domain.CheckIsometryofInvVectorTransport(&InitialX);
	//Domain.CheckVecTranComposeInverseVecTran(&InitialX);
	//Domain.CheckTranHInvTran(&InitialX);
	//return;

	//InitialX.Print("initialX:");

	// Generate the matrices in the Low rank approximation problem.
	double *y = new double[L * 2 + K * L * 2 + L * N * 2];
	double *B = y + L * 2;
	double *C = B + K * L * 2;
	for (integer i = 0; i < L * 2 + K * L * 2 + L * N * 2; i++)
		y[i] = genrandnormal();
	integer nzmaxB = L * K;
	integer nzmaxC = L * N;
	int *irB = new int[2 * L * K + 2 * L * N];
	int *jcB = irB + L * K;
	int *irC = jcB + L * K;
	int *jcC = irC + L * N;
	for (integer i = 0; i < K; i++)
	{
		for (integer j = 0; j < L; j++)
		{
			irB[j + i * L] = j;
			jcB[j + i * L] = i;
		}
	}
	for (integer i = 0; i < N; i++)
	{
		for (integer j = 0; j < L; j++)
		{
			irC[j + i * L] = j;
			jcC[j + i * L] = i;
		}
	}

	//blas_sparse_matrix sB;
	//sB = BLAS_duscr_begin(2 * L, K);
	////BLAS_duscr_insert_entries(sB, nzmaxB, B, irB, jcB);
	////BLAS_duscr_end(sB);
	//BLAS_usds(sB);

	CFR2BlindDeconvolution Prob(y, B, nzmaxB, irB, jcB, true, C, nzmaxC, irC, jcC, true, L, K, N, r, 0, 1, 1);
	Prob.SetDomain(&Domain);

	//Prob.CheckGradHessian(&InitialX);

	LRBFGS *RSDsolver = new LRBFGS(&Prob, &InitialX);
	//RSD *RSDsolver = new RSD(&Prob, &InitialX);
	//->LineSearch_LS = ARMIJO;
	//RSDsolver->LS_beta = 0.01;
	//RSDsolver->RCGmethod = DAI_YUAN;
	RSDsolver->Debug = FINALRESULT;
	RSDsolver->OutputGap = 1;
	RSDsolver->Max_Iteration = 100;
	RSDsolver->Accuracy = 1e-6;
	RSDsolver->Finalstepsize = 1;
	RSDsolver->Tolerance = 1e-6;
	RSDsolver->LengthSY = 0;
	RSDsolver->nu = 0;
	RSDsolver->LS_ratio1 = 0.3;
	RSDsolver->LS_ratio2 = 0.3;
	RSDsolver->InitSteptype = ONESTEP;
	//RSDsolver->CheckParams();
	RSDsolver->Run();
	//Prob.CheckGradHessian(&InitialX);//--
	//Prob.CheckGradHessian(RSDsolver->GetXopt());//--

	if (RSDsolver->Getnormgfgf0() < 1e-6)
		printf("SUCCESS!\n");
	else
		printf("FAIL!\n");
	delete RSDsolver;

	delete[] y;
	delete[] irB;

};

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 10)
	{
		mexErrMsgTxt("The number of arguments should be at least ten.\n");
	}
	double *yr, *yi, *y, *Br, *Bi, *B, *Cr, *Ci, *C, *X, *Xopt;
	yr = mxGetPr(prhs[0]);
	yi = mxGetPi(prhs[0]);
	Br = mxGetPr(prhs[1]);
	Bi = mxGetPi(prhs[1]);
	Cr = mxGetPr(prhs[2]);;
	Ci = mxGetPi(prhs[2]);
	X = mxGetPr(prhs[3]);
	/* dimensions of input matrices */
	integer L, K, N, HasHHR, r;
	L = mxGetM(prhs[0]);
	K = mxGetN(prhs[1]);
	N = mxGetN(prhs[2]);
	bool isWaveLet = (N == 1);
	if (mxGetM(prhs[1]) != L || (mxGetM(prhs[2]) != L && !isWaveLet))
	{
		mexErrMsgTxt("The size of B or the size of C is not correct.\n");
	}
	if (isWaveLet)
	{
		N = static_cast<integer> (mxGetScalar(prhs[2]));
		C = nullptr;
	}
	r = static_cast<integer> (mxGetScalar(prhs[4]));
	double rho, d, mu;
	rho = mxGetScalar(prhs[7]);
	d = mxGetScalar(prhs[8]);
	mu = mxGetScalar(prhs[9]);
	HasHHR = static_cast<integer> (mxGetScalar(prhs[5]));
	if (mxGetM(prhs[3]) != 2 * (K + N) * r)
	{
		mexErrMsgTxt("The size of initial x is not correct.\n");
	}

	bool isBsparse = mxIsSparse(prhs[1]);
	bool isCsparse = mxIsSparse(prhs[2]);
	integer nzmaxB = 0, nzmaxC = 0;
	size_t *irB = nullptr, *jcB = nullptr, *irC = nullptr, *jcC = nullptr;
	int *inirB = nullptr, *injcB = nullptr, *inirC = nullptr, *injcC = nullptr;
	if (isBsparse)
	{
		nzmaxB = mxGetNzmax(prhs[1]);
		irB = mxGetIr(prhs[1]);
		jcB = mxGetJc(prhs[1]);
		B = new double[2 * nzmaxB];

		for (integer i = 0; i < nzmaxB; i++)
		{
			B[2 * i] = Br[i];
			if (Bi == nullptr)
				B[2 * i + 1] = 0;
			else
				B[2 * i + 1] = Bi[i];
		}

		inirB = new int[2 * nzmaxB];
		injcB = inirB + nzmaxB;
		for (integer i = 0; i < K; i++)
		{
			for (unsigned long long j = jcB[i]; j < jcB[i + 1]; j++)
			{
				/*row: ir[j], column: i, entry: A[j]*/
				inirB[j] = irB[j];
				injcB[j] = i;
			}
		}
	}
	else
	{
		B = new double[2 * L * K];
		for (integer i = 0; i < L * K; i++)
		{
			B[2 * i] = Br[i];
			if (Bi == nullptr)
				B[2 * i + 1] = 0;
			else
				B[2 * i + 1] = Bi[i];
		}
	}
	if (isCsparse && !isWaveLet)
	{
		nzmaxC = mxGetNzmax(prhs[2]);
		irC = mxGetIr(prhs[2]);
		jcC = mxGetJc(prhs[2]);
		C = new double[2 * nzmaxC];

		for (integer i = 0; i < nzmaxC; i++)
		{
			C[2 * i] = Cr[i];
			if (Ci == nullptr)
				C[2 * i + 1] = 0;
			else
				C[2 * i + 1] = Ci[i];
		}

		inirC = new int[2 * nzmaxC];
		injcC = inirC + nzmaxC;
		for (integer i = 0; i < N; i++)
		{
			for (unsigned long long j = jcC[i]; j < jcC[i + 1]; j++)
			{
				/*row: ir[j], column: i, entry: A[j]*/
				inirC[j] = irC[j];
				injcC[j] = i;
			}
		}
	}
	else
	if(!isWaveLet)
	{
		C = new double[2 * L * N];
		for (integer i = 0; i < L * N; i++)
		{
			C[2 * i] = Cr[i];
			if (Ci == nullptr)
				C[2 * i + 1] = 0;
			else
				C[2 * i + 1] = Ci[i];
		}
	}

	y = new double[2 * L];
	for (integer i = 0; i < L; i++)
	{
		y[2 * i] = yr[i];
		if (yi == nullptr)
			y[2 * i + 1] = 0;
		else
			y[2 * i + 1] = yi[i];
	}

	genrandseed(0);
	CheckMemoryDeleted = new std::map<integer *, integer>;

	// Obtain an initial iterate
	CFR2Variable LRX(K, N, r);
	double *LRXptr = LRX.ObtainWriteEntireData();
	for (integer i = 0; i < 2 * (K + N) * r; i++)
		LRXptr[i] = X[i];

	CFR2Variable *LRsoln = nullptr;
	if (nrhs >= 8)
	{
		double *soln = mxGetPr(prhs[7]); /*soln: n by p*/
		LRsoln = new CFR2Variable(K, N, r);
		double *LRsolnptr = LRsoln->ObtainWriteEntireData();
		for (integer i = 0; i < 2 * (K + N) * r; i++)
		{
			LRsolnptr[i] = soln[i];
		}
	}

	// Define the manifold
	CFixedRank2Factors Domain(K, N, r);

	// Define the matrix completion problem
	CFR2BlindDeconvolution Prob(y, B, nzmaxB, inirB, injcB, isBsparse, C, nzmaxC, inirC, injcC, isCsparse, L, K, N, r, rho, d, mu);
	Prob.SetDomain(&Domain);

	Domain.SetHasHHR(HasHHR != 0);
	//Domain.CheckParams();

	// Call the function defined in DriverMexProb.h
	ParseSolverParamsAndOptimizing(prhs[6], &Prob, &LRX, LRsoln, plhs);

	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	delete LRsoln;
	delete[] inirB;
	delete[] inirC;
	delete[] y;
	delete[] B;
	delete[] C;

	return;
}

#endif
#endif
