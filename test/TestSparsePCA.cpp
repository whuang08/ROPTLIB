
#include "test/TestSparsePCA.h"

using namespace ROPTLIB;

#if !defined(MATLAB_MEX_FILE) && defined(TESTSPARSEPCA)

std::map<integer *, integer> *CheckMemoryDeleted;

int main(void)
{
	genrandseed(0);
	CheckMemoryDeleted = new std::map<integer *, integer>;

	integer p = 10, r = 3, n = 5;
	double mu = 5e-2;

	double *B = new double[p * n + r + 2 * p * r];
	double *Dsq = B + n * p;
	double *X = Dsq + r;
	double *Xopt = X + p * r;
	for (integer i = 0; i < p * n; i++)
		B[i] = genrandnormal();

	integer minpn = (p > n) ? n : p;
	integer maxpn = (p > n) ? p : n;
	double *Bcopy = new double[p * n + minpn + p * minpn + minpn * n];
	double *S = Bcopy + p * n;
	double *U = S + minpn;
	double *Vt = U + p * minpn;

	integer *iwork = new integer[8 * minpn];
	integer lwork = -1;
	double workoptsize = 0;
	integer np = n * p, inc = 1, unref = 1, info;
	// Bcopy <- B, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
	dcopy_(&np, B, &inc, Bcopy, &inc);
	char *jobs = const_cast<char *> ("s");
	// obtain optimal size for work
	printf("start obtaining optimal size for dgesdd\n");
	dgesdd_(jobs, &p, &n, Bcopy, &p, S, U, &p, Vt, &minpn, &workoptsize, &lwork, iwork, &info);
	lwork = static_cast<integer> (workoptsize);
	printf("the optimal size is %d\n", lwork);
	double *work = new double[lwork];
	// compute SVD for matrix Bcopy, such that Bcopy = U * diag(S) * Vt,
	// details: http://www.netlib.org/lapack/explore-html/db/db4/dgesdd_8f.html
	dgesdd_(jobs, &p, &n, Bcopy, &p, S, U, &p, Vt, &minpn, work, &lwork, iwork, &info);
	delete[] iwork;
	delete[] work;
	printf("info of dgesdd:%d\n", info);
	if (info != 0)
	{ // dgesdd fails, use dgesvd instead
		// Bcopy <- B, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		dcopy_(&np, B, &inc, Bcopy, &inc);
		lwork = -1;
		printf("start obtaining optimal size for dgesvd\n");
		dgesvd_(jobs, jobs, &p, &n, Bcopy, &p, S, U, &p, Vt, &minpn, &workoptsize, &lwork, &info);
		lwork = static_cast<integer> (workoptsize);
		printf("the optimal size is %d\n", lwork);
		work = new double[lwork];
		// compute SVD for matrix Bcopy, such that Bcopy = U * diag(S) * Vt,
		// details: http://www.netlib.org/lapack/explore-html/d8/d2d/dgesvd_8f.html
		dgesvd_(jobs, jobs, &p, &n, Bcopy, &p, S, U, &p, Vt, &minpn, work, &lwork, &info);
		printf("info of dgesvd:%d\n", info);
		delete[] work;
		if (info != 0)
		{
			printf("dgesvd fails!\n");
			return 0;
		}
	}
	for (integer i = 0; i < r; i++)
		Dsq[i] = S[i] * S[i];
	for (integer i = 0; i < p * r; i++)
	{
		X[i] = U[i];
	}
	delete[] Bcopy;

	testSparsePCA(B, Dsq, p, n, r, mu, X, Xopt);
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	delete[] B;

#ifdef _WIN64
#ifdef _DEBUG
	_CrtDumpMemoryLeaks();
#endif
#endif
	return 0;
}
#endif

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs < 3)
    {
        mexErrMsgTxt("The number of arguments should be at least three.\n");
    }
	double *B, *D, *X, *Xopt;
	B = mxGetPr(prhs[0]);
	D = mxGetPr(prhs[1]);
	X = mxGetPr(prhs[2]);
	double mu = mxGetScalar(prhs[3]);
	/* dimensions of input matrices */
	integer p, n, r;
	p = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);
	r = mxGetM(prhs[1]);
    if(mxGetN(prhs[1]) != 1)
    {
        mexErrMsgTxt("The size of matrix Dsq is not correct.\n");
    }
    if(mxGetM(prhs[2]) != p || mxGetN(prhs[2]) != r)
    {
        mexErrMsgTxt("The size of matrix X is not correct.\n");
    }
    
	printf("(p, n, r):%d,%d,%d\n", p, n, r);

	/*create output matrix*/
	plhs[0] = mxCreateDoubleMatrix(p, r, mxREAL);
	Xopt = mxGetPr(plhs[0]);

	genrandseed(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;
	testSparsePCA(B, D, p, n, r, mu, X, Xopt);
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

void testSparsePCA(double *B, double *Dsq, integer p, integer n, integer r, double mu, double *X, double *Xopt)
{
	// reduce n to r
	ObliqueVariable ObliqueX(p, r);

	Oblique Mani(p, r);
	Mani.CheckParams();
	//Mani.CheckIntrExtr(&ObliqueX);
	//Mani.CheckRetraction(&ObliqueX);
	//Mani.CheckDiffRetraction(&ObliqueX);
	//Mani.CheckLockingCondition(&ObliqueX);
	//Mani.CheckcoTangentVector(&ObliqueX);
	//Mani.CheckIsometryofVectorTransport(&ObliqueX);
	//Mani.CheckIsometryofInvVectorTransport(&ObliqueX);
	//Mani.CheckVecTranComposeInverseVecTran(&ObliqueX);
	//Mani.CheckTranHInvTran(&ObliqueX);
	//Mani.CheckHaddScaledRank1OPE(&ObliqueX);

	double *ObliqueXptr = ObliqueX.ObtainWriteEntireData();
	for (integer i = 0; i < ObliqueX.Getlength(); i++)
		ObliqueXptr[i] = X[i];
	ObliqueX.Print("initialX:");//---
	ObliqueSparsePCA Prob(B, Dsq, mu, p, n, r);
	Prob.SetDomain(&Mani);

	//Prob.CheckGradHessian(&ObliqueX);

	// test LRBFGS
	printf("********************************Check all line search algorithm in LRBFGSLPsub*************************************\n");
	RBFGSLPSub *LRBFGSLPSubsolver = new RBFGSLPSub(&Prob, &ObliqueX);
	LRBFGSLPSubsolver->Debug = ITERRESULT;
	//LRBFGSLPSubsolver->Tolerance = 1e-12;
	//LRBFGSLPSubsolver->InitSteptype = QUADINTMOD;
	LRBFGSLPSubsolver->NumExtraGF = Mani.GetIntrDim() * 5;
	LRBFGSLPSubsolver->lambdaLower = 1e-3;
	LRBFGSLPSubsolver->lambdaUpper = 1e3;
	LRBFGSLPSubsolver->Eps = 1e-3;
	LRBFGSLPSubsolver->Max_Iteration = 1000;
	LRBFGSLPSubsolver->OutputGap = 10;
	//LRBFGSLPSubsolver->CheckParams();
	LRBFGSLPSubsolver->Run();
	const Element *xopt = LRBFGSLPSubsolver->GetXopt();
	printf("[0:0.0005)  :%d\n", GetNumberBetweenC1andC2(xopt, 0, 0.0005));
	printf("[0.0005:0.1):%d\n", GetNumberBetweenC1andC2(xopt, 0.0005, 0.1));
	printf("[0.1:1)    :%d\n", GetNumberBetweenC1andC2(xopt, 0.1, 1));
	xopt->Print("Xopt:");//---
	if (Xopt != nullptr)
	{
		const double *xoptptr = xopt->ObtainReadData();
		for (integer i = 0; i < xopt->Getlength(); i++)
			Xopt[i] = xoptptr[i];
	}
	delete LRBFGSLPSubsolver;
};

integer GetNumberBetweenC1andC2(const Element *x, double c1, double c2)
{
	const double *xptr = x->ObtainReadData();
	integer length = x->Getlength();
	integer result = 0;
	for (integer i = 0; i < length; i++)
	{
		if (fabs(xptr[i]) >= c1 && fabs(xptr[i]) < c2)
		{
			result++;
		}
	}
	return result;
};
