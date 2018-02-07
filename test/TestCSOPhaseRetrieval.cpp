#include "test/TestCSOPhaseRetrieval.h"

#ifdef ROPTLIB_WITH_FFTW

using namespace ROPTLIB;

void testCSOPhaseRetrieval(void)
{
	integer n1 = 64, n2 = 64, r = 1, l = 6;
	integer n = n1 * n2, m = n * l;
	double kappa = 0.2;

	// Generate the matrices in the Low rank approximation problem.
	double *b = new double[3 * m];
	double *masks = b + m;
	for (integer i = 0; i < 3 * m; i++)
		b[i] = genrandreal();
	for (integer i = 0; i < 2 * m; i++)
		masks[i] = genrandnormal();

	double *xtrue = new double[2 * n];
	for (integer i = 0; i < 2 * n; i++)
		xtrue[i] = genrandnormal();

	double *ZY = new double[2 * m * r];
	double sqn = sqrt(static_cast<double> (n1 * n2));
	double *sZqi = new double[n1 * n2];

	for (integer i = 0; i < l; i++)
	{
		for (integer j = 0; j < n1 * n2; j++)
			sZqi[j] = 0;

		for (integer k = 0; k < n1 * n2; k++)
		{
			ZY[i * 2 * n1 * n2 + 2 * k] = (xtrue[2 * k] * masks[2 * k + i * 2 * n1 * n2] - xtrue[2 * k + 1] * masks[2 * k + 1 + i * 2 * n1 * n2]) / sqn;
			ZY[i * 2 * n1 * n2 + 2 * k + 1] = (xtrue[2 * k + 1] * masks[2 * k + i * 2 * n1 * n2] + xtrue[2 * k] * masks[2 * k + 1 + i * 2 * n1 * n2]) / sqn;
		}

		fftw_plan p = fftw_plan_dft_2d(n1, n2, (fftw_complex *)(ZY + i * 2 * n1 * n2), (fftw_complex *)(ZY + i * 2 * n1 * n2), FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(p);
		for (integer k = 0; k < n1 * n2; k++)
		{
			sZqi[k] = sZqi[k] + ZY[i * 2 * n1 * n2 + 2 * k] * ZY[i * 2 * n1 * n2 + 2 * k] + ZY[i * 2 * n1 * n2 + 2 * k + 1] * ZY[i * 2 * n1 * n2 + 2 * k + 1];
		}

		for (integer j = 0; j < n1 * n2; j++)
		{
			b[j + i * n1 * n2] = sZqi[j];
		}
	}
	delete[] sZqi;
	delete[] ZY;
	delete[] xtrue;

	double *soln = new double[2 * n1 * n2 + 2 * n1 * n2 * r];
	double *initX = soln + 2 * n1 * n2;
	for (integer i = 0; i < 2 * n1 * n2 * r; i++)
		initX[i] = genrandnormal();
	double time;
	integer nfft, nn;

	CSOPRRankReduce(initX, b, masks, n1, n2, l, r, kappa, 1000, soln, time, nfft, nn);

	//WFlow(initX, b, masks, n1, n2, l, r, 1000, soln, time, nfft);

	delete[] soln;
	delete[] b;
}

double delta = 0.95;
double tol = 1e-6;

void CSOPRRankReduce(double *initX, double *b, double *masks, integer n1, integer n2, integer l, integer r, double kappa, integer maxiter, double *outsoln, double &outtime, integer &outnfft, integer &outnn)
{
	integer n = n1 * n2, m = n * l;
	outtime = 0;
	outnfft = 0;
	outnn = 0;

	CSOPhaseRetrieval *Prob = nullptr;
	LRBFGS *solver = nullptr;
	CpxNStQOrth *Domain = nullptr;
	CSOVariable *InitialX = nullptr;
	InitialX = new CSOVariable(n, r);
	double *initXptr = InitialX->ObtainWriteEntireData();
	integer length = 2 * n * r;
	dcopy_(&length, initX, &GLOBAL::IONE, initXptr, &GLOBAL::IONE);

	double *S = new double[r + 2 * n * r];
	double *U = S + r;
	while (true)
	{
		printf("Rank is %d\n", r);
		kappa = (r == 1) ? 0 : kappa;

		Domain = new CpxNStQOrth(n, r);
		Domain->SetHasHHR(false);
		Prob = new CSOPhaseRetrieval(b, masks, kappa, n1, n2, l, r);
		Prob->SetDomain(Domain);

		solver = new LRBFGS(Prob, InitialX);
		solver->Debug = FINALRESULT;
		solver->OutputGap = 5;
		solver->Max_Iteration = maxiter;
		solver->Min_Iteration = 10;
		solver->StopPtr = &MyPRstop;
		solver->BBratio = 0;
		//solver->CheckParams();
		//solver->LengthSY = 2;
		solver->Accuracy = 1e-6;
		solver->Finalstepsize = 1;
		solver->Tolerance = 1e-10;
		solver->Run();

		Variable *Xopt = solver->GetXopt()->ConstructEmpty();
		solver->GetXopt()->CopyTo(Xopt);
		outtime += solver->GetComTime();
		outnfft += (solver->Getnf() + solver->Getng()) * l * r;
		outnn += 2 * r * r * solver->Getng() + 6 * r * r * solver->GetIter() + ((r > 1) ? 2 * (solver->GetIter() - 1) * r * r : 0);

		double *solnptr = Xopt->ObtainWritePartialData();
		delete solver;
		delete Prob;
		delete InitialX;
		InitialX = nullptr;
		delete Domain;

		if (r > 1)
		{
			doublecomplex workoptsize;
			integer lwork = -1, info;
			double *rwork = new double[10 * r];
#ifndef MATLAB_MEX_FILE
			zgesvd_(GLOBAL::S, GLOBAL::N, &n, &r, (doublecomplex *)solnptr, &n, S, (doublecomplex *)U, &n, (doublecomplex *)U, &n, &workoptsize, &lwork, rwork, &info);
#else
			zgesvd_(GLOBAL::S, GLOBAL::N, &n, &r, const_cast<double *> (solnptr), &n, S, U, &n, U, &n, (double *) &workoptsize, &lwork, rwork, &info);
#endif
			lwork = static_cast<integer> (workoptsize.r);
			double *work = new double[2 * lwork];
#ifndef MATLAB_MEX_FILE
			zgesvd_(GLOBAL::S, GLOBAL::N, &n, &r, (doublecomplex *)solnptr, &n, S, (doublecomplex *)U, &n, (doublecomplex *)U, &n, (doublecomplex *) work, &lwork, rwork, &info);
#else
			zgesvd_(GLOBAL::S, GLOBAL::N, &n, &r, const_cast<double *> (solnptr), &n, S, U, &n, U, &n, work, &lwork, rwork, &info);
#endif
			if (info != 0)
			{
				printf("Error:singular value decomposition failed!\n");
			}
			delete[] work;
			delete[] rwork;
			double normSdsqrk = sqrt(ddot_(&r, S, &GLOBAL::IONE, S, &GLOBAL::IONE) / r);
			for (integer i = 0; i < r; i++)
			{
				if (S[i] / normSdsqrk < delta)
				{
					InitialX = new CSOVariable(n, i);
					double *InitXptr = InitialX->ObtainWriteEntireData();
					integer length = 2 * n * i, n2 = n * 2;
					dcopy_(&length, const_cast<double *> (solnptr), &GLOBAL::IONE, InitXptr, &GLOBAL::IONE);
					for(integer j = 0; j < i; j++)
						dscal_(&n2, S + j, InitXptr + n2 * j, &GLOBAL::IONE);
					break;
				}
			}
		}

		if (r == 1)
		{
			/*find rank-1 solution*/
			for (integer i = 0; i < 2 * n; i++)
			{
				outsoln[i] = solnptr[i];
			}
			delete Xopt;
			break;
		}
		delete Xopt;
		if (InitialX == nullptr)
		{
			r--;
			InitialX = new CSOVariable(n, r);
			InitialX->RandInManifold();
		}
		else
		{
			r = InitialX->Getsize()[1];
		}
	}
	delete[] S;
	/*finaliterate, time, nfft, nn*/
};

bool MyPRstop(Variable *x, Vector *gf, double f, double ngf, double ngf0, const Problem *prob, const Solvers *solver)
{
	integer r = x->Getsize()[1];
	bool result = false;
	if (r > 1)
	{
		integer n = x->Getsize()[0] / 2;
		Variable *tmp = x->ConstructEmpty();
		x->CopyTo(tmp);
		double *solnptr = tmp->ObtainWritePartialData();

		double *S = new double[r];
		doublecomplex workoptsize, uv;
		integer lwork = -1, info;
		double *rwork = new double[10 * r];
#ifndef MATLAB_MEX_FILE
		zgesvd_(GLOBAL::N, GLOBAL::N, &n, &r, (doublecomplex *)solnptr, &n, S, &uv, &GLOBAL::IONE, &uv, &GLOBAL::IONE, &workoptsize, &lwork, rwork, &info);
#else
		zgesvd_(GLOBAL::N, GLOBAL::N, &n, &r, solnptr, &n, S, (double *) &uv, &GLOBAL::IONE, (double *) &uv, &GLOBAL::IONE, (double *) &workoptsize, &lwork, rwork, &info);
#endif
		lwork = static_cast<integer> (workoptsize.r);
		double *work = new double[2 * lwork];
#ifndef MATLAB_MEX_FILE
		zgesvd_(GLOBAL::N, GLOBAL::N, &n, &r, (doublecomplex *)solnptr, &n, S, &uv, &GLOBAL::IONE, &uv, &GLOBAL::IONE, (doublecomplex *)work, &lwork, rwork, &info);
#else
		zgesvd_(GLOBAL::N, GLOBAL::N, &n, &r, solnptr, &n, S, (double *) &uv, &GLOBAL::IONE, (double *)&uv, &GLOBAL::IONE, work, &lwork, rwork, &info);
#endif
		if (info != 0)
		{
			printf("Error:singular value decomposition failed!\n");
		}
		delete[] work;
		delete[] rwork;
		delete tmp;
		double normSdsqrk = sqrt(ddot_(&r, S, &GLOBAL::IONE, S, &GLOBAL::IONE) / r);
		for (integer i = 0; i < r; i++)
		{
			if (S[i] / normSdsqrk < delta)
			{
				result = true;
				break;
			}
		}
		delete[] S;
		return result;
	}

	return ngf < tol;
};

void WFegf(double *x, double *b, double *masks, double *egf, integer n1, integer n2, integer l, integer r)
{
	integer m = n1 * n2 * l;
	double *cD = new double[m];
	double *ZY = new double[2 * m * r];
	double sqn = sqrt(static_cast<double> (n1 * n2));

	for (integer i = 0; i < m; i++)
		cD[i] = 0;

	double *sZqi = new double[n1 * n2];

	for (integer i = 0; i < l; i++)
	{
		for (integer j = 0; j < n1 * n2; j++)
			sZqi[j] = 0;

		for (integer j = 0; j < r; j++)
		{
			for (integer k = 0; k < n1 * n2; k++)
			{
				ZY[i * 2 * n1 * n2 + 2 * k + j * 2 * m] = (x[2 * k + j * 2 * n1 * n2] * masks[2 * k + i * 2 * n1 * n2] - x[2 * k + 1 + j * 2 * n1 * n2] * masks[2 * k + 1 + i * 2 * n1 * n2]) / sqn;
				ZY[i * 2 * n1 * n2 + 2 * k + 1 + j * 2 * m] = (x[2 * k + 1 + j * 2 * n1 * n2] * masks[2 * k + i * 2 * n1 * n2] + x[2 * k + j * 2 * n1 * n2] * masks[2 * k + 1 + i * 2 * n1 * n2]) / sqn;
			}

			fftw_plan p = fftw_plan_dft_2d(n1, n2, (fftw_complex *)(ZY + i * 2 * n1 * n2 + j * 2 * m), (fftw_complex *)(ZY + i * 2 * n1 * n2 + j * 2 * m), FFTW_FORWARD, FFTW_ESTIMATE);
			fftw_execute(p);

			for (integer k = 0; k < n1 * n2; k++)
			{
				sZqi[k] = sZqi[k] + ZY[i * 2 * n1 * n2 + 2 * k + j * 2 * m] * ZY[i * 2 * n1 * n2 + 2 * k + j * 2 * m] + ZY[i * 2 * n1 * n2 + 2 * k + 1 + j * 2 * m] * ZY[i * 2 * n1 * n2 + 2 * k + 1 + j * 2 * m];
			}
		}

		for (integer j = 0; j < n1 * n2; j++)
		{
			cD[j + i * n1 * n2] = sZqi[j] - b[j + i * n1 * n2];
		}
	}
	delete[] sZqi;

	double *DZY = new double[2 * m * r + 2 * n1 * n2 * r];
	double *tmp = DZY + 2 * m * r;
	for (integer i = 0; i < r; i++)
	{
		for (integer j = 0; j < m; j++)
		{
			DZY[2 * j + i * 2 * m] = cD[j] * ZY[2 * j + i * 2 * m];
			DZY[2 * j + 1 + i * 2 * m] = cD[j] * ZY[2 * j + 1 + i * 2 * m];
		}
	}
	delete[] cD;
	delete[] ZY;

	for (integer i = 0; i < 2 * n1 * n2 * r; i++)
		egf[i] = 0;

	integer length = 2 * n1 * n2 * r;
	for (integer i = 0; i < l; i++)
	{
		for (integer j = 0; j < r; j++)
		{
			fftw_plan p = fftw_plan_dft_2d(n1, n2, (fftw_complex *)(DZY + i * 2 * n1 * n2 + j * 2 * m), (fftw_complex *)(DZY + i * 2 * n1 * n2 + j * 2 * m), FFTW_BACKWARD, FFTW_ESTIMATE);
			fftw_execute(p);
			
			for (integer k = 0; k < n1 * n2; k++)
			{
				double realv = (DZY[2 * k + i * 2 * n1 * n2 + j * 2 * m] * masks[2 * k + i * 2 * n1 * n2] + DZY[2 * k + 1 + i * 2 * n1 * n2 + j * 2 * m] * masks[2 * k + 1 + i * 2 * n1 * n2]) / sqn;
				double imagv = (-DZY[2 * k + i * 2 * n1 * n2 + j * 2 * m] * masks[2 * k + 1 + i * 2 * n1 * n2] + DZY[2 * k + 1 + i * 2 * n1 * n2 + j * 2 * m] * masks[2 * k + i * 2 * n1 * n2]) / sqn;
				tmp[2 * k + j * n1 * n2 * 2] = realv;
				tmp[2 * k + 1 + j * n1 * n2 * 2] = imagv;
			}
		}
		daxpy_(&length, &GLOBAL::DONE, tmp, &GLOBAL::IONE, egf, &GLOBAL::IONE);
	}
	delete[] DZY;

	double scalar = 4.0 / ddot_(&m, b, &GLOBAL::IONE, b, &GLOBAL::IONE);
	dscal_(&length, &scalar, egf, &GLOBAL::IONE);
};

void WFlow(double *initX, double *b, double *masks, integer n1, integer n2, integer l, integer r, integer maxiter, double *outsoln, double &outtime, integer &outnfft)
{
	integer n = n1 * n2, m = n * l;
	double *egf = new double[2 * n * r];
	WFegf(initX, b, masks, egf, n1, n2, l, r);
	integer length = 2 * n * r;
	double ngf = sqrt(ddot_(&length, egf, &GLOBAL::IONE, egf, &GLOBAL::IONE)), ngf0 = ngf;
	
	unsigned long starttime = getTickCount();
	integer iter = 0;
	double deno = ddot_(&length, initX, &GLOBAL::IONE, initX, &GLOBAL::IONE);
	double stepsize = 0.05;

	dcopy_(&length, initX, &GLOBAL::IONE, outsoln, &GLOBAL::IONE);
	while (ngf > tol && iter < maxiter)
	{
		/* this stepsize works well for n = 2 ^ 2 to 256 ^ 2*/
		if (l == 6)
		{
			if (n == 8 * 8)
			{
                stepsize = ((1.0 - std::exp(-(1.0 + iter) / 330.0) < 0.13) ? (1.0 - std::exp(-(1.0 + iter) / 330.0)) : 0.13) * 5 / deno;
			}
			else
			if (n == 16 * 16)
			{
                stepsize = ((1.0 - std::exp(-(1.0 + iter) / 330.0) < 0.15) ? (1.0 - std::exp(-(1.0 + iter) / 330.0)) : 0.15) * 20 / deno;
			}
			else
			if (n == 32 * 32)
			{
                stepsize = ((1.0 - std::exp(-(1.0 + iter) / 330.0) < 0.17) ? (1.0 - std::exp(-(1.0 + iter) / 330.0)) : 0.17) * 50 / deno;
			}
		}

		if (l == 20)
		{
            stepsize = ((1.0 - std::exp(-(1.0 + iter) / 330.0) < 0.036 * std::log(static_cast<double> (n)) - 0.015) ? (1.0 - std::exp(-(1.0 + iter) / 330.0)) : 0.036 * std::log(static_cast<double> (n)) - 0.015) * (std::exp(2.35 * pow(static_cast<double> (n), 0.125) - 1.7)) / deno;
		}
		WFegf(outsoln, b, masks, egf, n1, n2, l, r);
		ngf = sqrt(ddot_(&length, egf, &GLOBAL::IONE, egf, &GLOBAL::IONE));
		if ((iter % 100) == 0)
		{
			printf("iter:%d, ngf:%3.2e, ngf/ngf0:%3.2e, stepsize:%3.2e\n", iter, ngf, ngf/ngf0, stepsize);
		}
		double scalar = -stepsize;
		daxpy_(&length, &scalar, egf, &GLOBAL::IONE, outsoln, &GLOBAL::IONE);
		iter++;
	}
	printf("total iter:%d, ngf:%3.2e, ngf/ngf0:%3.2e, stepsize:%3.2e\n", iter, ngf, ngf / ngf0, stepsize);

	delete[] egf;

	outtime = static_cast<double>(getTickCount() - starttime) / CLK_PS;
	outnfft = (iter + 1) * 2 * l;

};

void testCSOPhaseRetrievalFixedRank(void)
{
	integer n1 = 4, n2 = 8, r = 2, l = 4;
	integer n = n1 * n2, m = n * l;
	double kappa = 0.2;

	CpxNStQOrth Domain(n, r);
	Domain.SetHasHHR(false);
	CSOVariable InitialX(n, r);
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
	double *b = new double[3 * m];
	double *masks = b + m;
	for (integer i = 0; i < 3 * m; i++)
		b[i] = genrandreal();
	for (integer i = 0; i < 2 * m; i++)
		masks[i] = genrandnormal();

	//ForDebug::Print("b:", b, m);//----
	//ForDebug::Print("masks:", masks, 2 * m);//---

	//InitialX.Print("InitialX:");///----

	CSOPhaseRetrieval Prob(b, masks, kappa, n1, n2, l, r);
	Prob.SetDomain(&Domain);

	//std::cout << "f:" << Prob.f(&InitialX) << std::endl;//----

	//Vector *gf = Domain.GetEMPTYINTR()->ConstructEmpty();
	//Prob.Grad(&InitialX, gf);
	//delete gf;

	//Prob.CheckGradHessian(&InitialX);

	LRBFGS *RSDsolver = new LRBFGS(&Prob, &InitialX);
	//->LineSearch_LS = ARMIJO;
	//RSDsolver->LS_beta = 0.01;
	//RSDsolver->RCGmethod = DAI_YUAN;
	RSDsolver->Debug = FINALRESULT;
	RSDsolver->OutputGap = 100;
	RSDsolver->Max_Iteration = 100;
	//RSDsolver->CheckParams();
	RSDsolver->Accuracy = 1e-6;
	RSDsolver->Finalstepsize = 1;
	RSDsolver->Tolerance = 1e-6;
	RSDsolver->Run();
	if (RSDsolver->Getnormgfgf0() < 1e-6)
		printf("SUCCESS!\n");
	else
		printf("FAIL!\n");
	//Prob.CheckGradHessian(&InitialX);//--
	//Prob.CheckGradHessian(RSDsolver->GetXopt());//--


	//// Compute the smallest eigenvalue of the Hessian at initial iterate.
	//SphereTx DomainPH0(Prob.GetDomain(), &InitialX);
	//SphereTxRQ ProbHess0(Prob.GetDomain(), &InitialX, &Prob, true);
	//ProbHess0.SetDomain(&DomainPH0);
	//Variable *TV00 = DomainPH0.RandominManifold();
	//RTRNewton *RTRNewtonsolver = new RTRNewton(&ProbHess0, TV00);
	//RTRNewtonsolver->Debug = NOOUTPUT;
	//RTRNewtonsolver->Run();
	//if (RTRNewtonsolver->Getnormgfgf0() > 1e-4)
	//	printf("Stop early when finding the smallest eigenvalue of the Hessian at initial iterate\n");

	//delete RTRNewtonsolver;
	//ProbHess0.SetMinorMax(false);
	//RTRNewtonsolver = new RTRNewton(&ProbHess0, TV00);
	//RTRNewtonsolver->Debug = NOOUTPUT;
	//RTRNewtonsolver->Run();
	//if (RTRNewtonsolver->Getnormgfgf0() > 1e-4)
	//	printf("Stop early when finding the largest eigenvalue of the Hessian at initial iterate\n");
	//delete TV00;

	//// Compute the smallest eigenvalue of the Hessian at root.
	//Variable *root = InitialX.ConstructEmpty();
	//RSDsolver->GetXopt()->CopyTo(root);
	//SphereTx DomainPH(Prob.GetDomain(), root);
	//SphereTxRQ ProbHess(Prob.GetDomain(), root, &Prob, true);
	//ProbHess.SetDomain(&DomainPH);
	//Variable *TV0 = DomainPH.RandominManifold();
	//delete RTRNewtonsolver;
	//RTRNewtonsolver = new RTRNewton(&ProbHess, TV0);
	//RTRNewtonsolver->Debug = NOOUTPUT;
	//RTRNewtonsolver->Run();
	//if (RTRNewtonsolver->Getnormgfgf0() > 1e-4)
	//	printf("Stop early when finding the smallest eigenvalue of the Hessian at optimum\n");

	//delete RTRNewtonsolver;
	//ProbHess.SetMinorMax(false);
	//RTRNewtonsolver = new RTRNewton(&ProbHess, TV0);
	//RTRNewtonsolver->Debug = NOOUTPUT;
	//RTRNewtonsolver->Run();
	//if (RTRNewtonsolver->Getnormgfgf0() > 1e-4)
	//	printf("Stop early when finding the largest eigenvalue of the Hessian at optimum\n");
	//delete RTRNewtonsolver;
	//delete root;
	//delete TV0;



	delete RSDsolver;

	delete[] b;
};

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 8)
	{
		mexErrMsgTxt("The number of arguments should be at least eight.\n");
	}
	double *b, *masks, *initX;
	bool isRie = true;
	integer maxiter = 0;
	double kappa;
	b = mxGetPr(prhs[0]);
	masks = mxGetPr(prhs[1]);
	initX = mxGetPr(prhs[2]);
	kappa = mxGetScalar(prhs[3]);
	delta = mxGetScalar(prhs[4]);
	tol = mxGetScalar(prhs[5]);
	isRie = (mxGetScalar(prhs[6]) != 0);
	maxiter = static_cast<integer> (mxGetScalar(prhs[7]));
	
	/* dimensions of input matrices */
	integer m, n1, n2, l, n, r;
	const size_t *size = mxGetDimensions(prhs[1]);
	n1 = size[0] / 2;
	n2 = size[1];
	l = size[2];
	m = mxGetM(prhs[0]);
	n = mxGetM(prhs[2]) / 2;
	r = mxGetN(prhs[2]);
	if (n != n1 * n2)
	{
		mexErrMsgTxt("The size of masks or the size of initX is not correct.\n");
	}
	if (m != n * l)
	{
		mexErrMsgTxt("The size of b or the size of masks is not correct.\n");
	}
	if (!isRie && r != 1)
	{
		mexErrMsgTxt("r need be 1 for WF method.\n");
	}

	genrandseed(0);
	CheckMemoryDeleted = new std::map<integer *, integer>;

	plhs[0] = mxCreateDoubleMatrix(2 * n1 * n2, 1, mxREAL);

	double *soln = mxGetPr(plhs[0]);
	double time = 0;
	integer nfft = 0, nn = 0;
	if(isRie)
		CSOPRRankReduce(initX, b, masks, n1, n2, l, r, kappa, maxiter, soln, time, nfft, nn);
	else
		WFlow(initX, b, masks, n1, n2, l, r, maxiter, soln, time, nfft);

	plhs[1] = mxCreateDoubleScalar(static_cast<double> (time));
	plhs[2] = mxCreateDoubleScalar(static_cast<double> (nfft));
	plhs[3] = mxCreateDoubleScalar(static_cast<double> (nn));

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
#endif

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
