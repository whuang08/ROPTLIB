#include "test/TestCFR2BlindDecon2D.h"

#ifdef ROPTLIB_WITH_FFTW

using namespace ROPTLIB;


void testCFR2BlindDecon2D(void)
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

	//integer nn1 = 4, nn2 = 4;
	//doublecomplex *vv = new doublecomplex[4 * 4];
	//for (integer i = 0; i < nn1 * nn2; i++)
	//{
	//	vv[i].r = i;
	//	vv[i].i = i;
	//}
	//std::cout << "original:" << std::endl;//---
	//for (integer i = 0; i < nn1; i++)
	//{
	//	for (integer j = 0; j < nn2; j++)
	//	{
	//		std::cout << vv[i + j * nn1].r << "+ i " << vv[i + j * nn1].i << "\t";
	//	}
	//	std::cout << std::endl;//---
	//}
	//haarFWT_2d(4, 4, vv);
	//std::cout << "Haar transform:" << std::endl;//---
	//for (integer i = 0; i < nn1; i++)
	//{
	//	for (integer j = 0; j < nn2; j++)
	//	{
	//		std::cout << vv[i + j * nn1].r << "+ i " << vv[i + j * nn1].i << "\t";
	//	}
	//	std::cout << std::endl;//---
	//}
	//std::cout << "Haar inverse transform:" << std::endl;//---
	//haarFWT_2d_inverse(nn1, nn2, vv);

	//for (integer i = 0; i < nn1; i++)
	//{
	//	for (integer j = 0; j < nn2; j++)
	//	{
	//		std::cout << vv[i + j * nn1].r << "+ i " << vv[i + j * nn1].i << "\t";
	//	}
	//	std::cout << std::endl;//---
	//}
	//return;

	integer n1 = 4, n2 = 4, r = 1, L = n1 * n2;

	CFixedRank2Factors Domain(L, L, r);
	Domain.SetHasHHR(false);
	CFR2Variable InitialX(L, L, r);
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
	double *y = new double[L * 2 + L * L * 2 + L * L * 2];
	double *B = y + L * 2;
	double *C = B + L * L * 2;
	for (integer i = 0; i < L * 2 + L * L * 2 + L * L * 2; i++)
		y[i] = genrandnormal();
	integer nzmaxB = L * L;
	integer nzmaxC = L * L;
	int *irB = new int[2 * L * L + 2 * L * L];
	int *jcB = irB + L * L;
	int *irC = jcB + L * L;
	int *jcC = irC + L * L;
	for (integer i = 0; i < L; i++)
	{
		for (integer j = 0; j < L; j++)
		{
			irB[j + i * L] = j;
			jcB[j + i * L] = i;
		}
	}
	for (integer i = 0; i < L; i++)
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
	CFR2BlindDecon2D Prob(y, B, nzmaxB, irB, jcB, true, C, nzmaxC, irC, jcC, true, n1, n2, r, 0, 1, 1);
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
	if (nrhs < 12)
	{
		mexErrMsgTxt("The number of arguments should be at least twelve.\n");
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
	integer L, K, N, HasHHR, n1, n2, r;
	L = mxGetM(prhs[0]);
	K = mxGetN(prhs[1]);
	N = mxGetN(prhs[2]);
	bool isWaveLet = (N == 1);
	if (K != L || mxGetM(prhs[1]) != L || (mxGetM(prhs[2]) != L && !isWaveLet))
	{
		mexErrMsgTxt("The size of B or the size of C is not correct.\n");
	}
	if (isWaveLet)
	{
		N = static_cast<integer> (mxGetScalar(prhs[2]));
		if(L != N)
		{
			mexErrMsgTxt("The size of C is not correct.\n");
		}

		C = nullptr;
	}
	n1 = static_cast<integer> (mxGetScalar(prhs[4]));
	n2 = static_cast<integer> (mxGetScalar(prhs[5]));
	r = static_cast<integer> (mxGetScalar(prhs[6]));
	double rho, d, mu;
	rho = mxGetScalar(prhs[9]);
	d = mxGetScalar(prhs[10]);
	mu = mxGetScalar(prhs[11]);
	HasHHR = static_cast<integer> (mxGetScalar(prhs[7]));
	if (mxGetM(prhs[3]) != 4 * L * r)
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
	if (nrhs >= 13)
	{
		double *soln = mxGetPr(prhs[12]); /*soln: n by p*/
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
	CFR2BlindDecon2D Prob(y, B, nzmaxB, inirB, injcB, isBsparse, C, nzmaxC, inirC, injcC, isCsparse, n1, n2, r, rho, d, mu);
	Prob.SetDomain(&Domain);

	Domain.SetHasHHR(HasHHR != 0);
	//Domain.CheckParams();

	// Call the function defined in DriverMexProb.h
	ParseSolverParamsAndOptimizing(prhs[8], &Prob, &LRX, LRsoln, plhs);

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
