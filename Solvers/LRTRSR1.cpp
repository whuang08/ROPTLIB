
#include "Solvers/LRTRSR1.h"

/*Define the namespace*/
namespace ROPTLIB{

	LRTRSR1::LRTRSR1(const Problem *prob, const Variable *initialx, const Variable *insoln)
	{
		Initialization(prob, initialx, insoln);
	};

	void LRTRSR1::SetProbX(const Problem *prob, const Variable *initialx, const Variable *insoln)
	{
		SolversTR::SetProbX(prob, initialx, insoln);
		const Vector *EMPTYETA;
		if (prob->GetDomain()->GetIsIntrinsic())
			EMPTYETA = prob->GetDomain()->GetEMPTYINTR();
		else
			EMPTYETA = prob->GetDomain()->GetEMPTYEXTR();
		s = EMPTYETA->ConstructEmpty();
		y = EMPTYETA->ConstructEmpty();

		prob->SetUseGrad(true);
		prob->SetUseHess(false);
	};

	void LRTRSR1::SetDefaultParams()
	{
		SolversTR::SetDefaultParams();
		theta = 0.1;
		kappa = 0.1;
		isconvex = false;
		LengthSY = 4;
		S = nullptr;
		Y = nullptr;
		YMGS = nullptr;
		inpss = 0;
		inpsy = 0;
		inpyy = 0;
		Currentlength = 0;
		beginidx = 0;
		SS = nullptr;
		SY = nullptr;
		PMGQ = nullptr;
		LU_PMGQ = nullptr;
		P = nullptr;
		gamma = 1;
		SolverName.assign("LRTRSR1");

		/*LSR solver tmp variables*/
		Psi = nullptr;
		Q = nullptr;
		RT = nullptr;
		U = nullptr;
		D = nullptr;
		P_parallel = nullptr;
		g_parallel = nullptr;
		Psig = nullptr;
		Lambda = nullptr;
		a_j = nullptr;
		tolLocalNewton = 1e-10;
	};

	LRTRSR1::~LRTRSR1(void)
	{
		delete s;
		delete y;
		DeleteVectors(S, LengthSY);
		DeleteVectors(Y, LengthSY);
		DeleteVectors(YMGS, LengthSY);
		if (SS != nullptr)
			delete[] SS;
		if (SY != nullptr)
			delete[] SY;
		if (PMGQ != nullptr)
			delete[] PMGQ;
		if (LU_PMGQ != nullptr)
			delete[] LU_PMGQ;
		if (P != nullptr)
			delete[] P;
		if (Psi != nullptr)
			delete[] Psi;
	};

	void LRTRSR1::Run(void)
	{
		DeleteVectors(S, LengthSY);
		NewVectors(S, LengthSY);
		DeleteVectors(Y, LengthSY);
		NewVectors(Y, LengthSY);
		DeleteVectors(YMGS, LengthSY);
		NewVectors(YMGS, LengthSY);
		if (SS != nullptr)
			delete[] SS;
		SS = new double[LengthSY * LengthSY];
		if (SY != nullptr)
			delete[] SY;
		SY = new double[LengthSY * LengthSY];
		if (PMGQ != nullptr)
			delete[] PMGQ;
		PMGQ = new double[LengthSY * LengthSY];
		if (LU_PMGQ != nullptr)
			delete[] LU_PMGQ;
		LU_PMGQ = new double[LengthSY * LengthSY];
		if (P != nullptr)
			delete[] P;
		P = new integer[LengthSY];

		if (Psi != nullptr)
			delete[] Psi;
		Psi = new double[3 * gf1->Getlength() * LengthSY + 2 * LengthSY * LengthSY + 5 * LengthSY + 2];
		Q = Psi + gf1->Getlength() * LengthSY;
		RT = Q + gf1->Getlength() * LengthSY;
		U = RT + LengthSY * LengthSY;
		D = U + LengthSY * LengthSY;
		P_parallel = D + LengthSY;
		g_parallel = P_parallel + gf1->Getlength() * LengthSY;
		Psig = g_parallel + LengthSY;
		Lambda = Psig + LengthSY;
		a_j = Lambda + LengthSY + 1;

		//double *Q = new double[2 * n * k + 2 * k * k + 3 * k];
		//double *RT = Q + n * k;
		//double *U = RT + k * k;
		//double *D = U + k * k;
		//double *P_parallel = D + k;
		//double *g_parallel = P_parallel + n * k;
		//double *Psig = g_parallel + k;
		//double *Lambda = new double[k + 1];
		//double *a_j = new double[k + 1];
		SolversTR::Run();
	};

	void LRTRSR1::tCG_TR(void)
	{/*This method is only used when the Euclidean metric or an intrinsic approach is used.*/

		integer n = gf1->Getlength(); /*number of rows*/
		integer k = Currentlength;    /*number of columns*/

		integer idx;
		/*Compute Psi = Y - gamma S, PMGS = SY - gamma SS */
		for (integer i = 0; i < Currentlength; i++)
		{
			idx = (i + beginidx) % LengthSY;
			Mani->scalarVectorAddVector(x1, -gamma, S[idx], Y[idx], YMGS[i]);
			dcopy_(&n, const_cast<double *> (YMGS[i]->ObtainReadData()), &GLOBAL::IONE, Psi + i * n, &GLOBAL::IONE);
		}
		for (integer i = 0; i < Currentlength; i++)
		{
			for (integer j = 0; j < Currentlength; j++)
			{
				PMGQ[i + j * Currentlength] = SY[i + j * LengthSY] - gamma * SS[i + j * LengthSY];
			}
		}
		/*
		QR decomposion:           Psi = Q RT^T
		eigenvalue decomposition: U D U^T = eig(RT^T * PMGQ^{-1} RT)
		P_parallel = Q * U
		Therefore: Psi * (PMGQ)^{-1} * Psi = Q * U * D * U^T * Q^T = P_parallel * D * P_parallel^T
		*/
		if (k > 0)
		{
			integer nk = n * k;
			dcopy_(&nk, Psi, &GLOBAL::IONE, Q, &GLOBAL::IONE);

			/*QR decomposition of Psi*/
			double *tau = new double[k];
			integer *jpvt = new integer[k];
			integer info;
			integer lwork = -1;
			double lworkopt;
			for (integer i = 0; i < k; i++)
				jpvt[i] = i + 1;
			// compute the space required in the dgeqp3
			dgeqp3_(&n, &k, Q, &n, jpvt, tau, &lworkopt, &lwork, &info);
			lwork = static_cast<integer> (lworkopt);
			double *work = new double[lwork];
			// QR decomposition for ptrHHR using Householder reflections. Householder reflectors and R are stored in ptrHHR.
			// details: http://www.netlib.org/lapack/explore-html/db/de5/dgeqp3_8f.html
			dgeqp3_(&n, &k, Q, &n, jpvt, tau, work, &lwork, &info);
			if (info < 0)
				printf("Error in qr decomposition!\n");
			for (integer i = 0; i < k; i++)
			{
				if (jpvt[i] != (i + 1))
					printf("Error in qf tCG in LRTRSR1!\n");
			}

			for (integer i = 0; i < k; i++)
			{
				for (integer j = 0; j < k; j++)
				{
					if (i > j)
						RT[j + i * k] = 0;
					else
						RT[j + i * k] = Q[i + j * n];
				}
			}
			// Generate an orthonormal matrix by using the Householder refections in resultM, output is stored in ptrHHR,
			// details: http://www.netlib.org/lapack/explore-html/d9/d1d/dorgqr_8f.html
			dorgqr_(&n, &k, &k, Q, &n, tau, work, &lwork, &info);
			if (info < 0)
				printf("Error in dorgqr of LRTRSR1::tCG_TR(void)!\n");

			delete[] work;
			delete[] jpvt;
			delete[] tau;

			double *tmp = new double[k * k];
			integer k2 = k * k;
			dcopy_(&k2, PMGQ, &GLOBAL::IONE, LU_PMGQ, &GLOBAL::IONE);
			dcopy_(&k2, RT, &GLOBAL::IONE, tmp, &GLOBAL::IONE);

			/* Compute U D U^T = eig(RT^T * PMGQ^{-1} RT) */
			// LU decomposion for PMGQ, PMGQ = P * L * U, L and U are stored in PMGQ, the permutation matrix is in P
			// details: http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html
			dgetrf_(&k, &k, LU_PMGQ, &k, P, &info);
			if (info < 0)
				printf("Error in forming Q matrix!\n");

			// solve linear system: PMGQ * X = v using the LU decomposition results from dgetrf, then solution is stored in v.
			// details: http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f.html
			dgetrs_(GLOBAL::N, &k, &k, LU_PMGQ, &k, P, tmp, &k, &info);

			dgemm_(GLOBAL::T, GLOBAL::N, &k, &k, &k, &GLOBAL::DONE, RT, &k, tmp, &k, &GLOBAL::DZERO, U, &k);
			//ForDebug::Print("RMR:", U, k, k);///---
			delete[] tmp;
			for (integer i = 0; i < k; i++)
			{
				for (integer j = i + 1; j < k; j++)
				{
					U[j + i * k] = (U[j + i * k] + U[i + j * k]) / 2;
					U[i + j * k] = U[j + i * k];
				}
			}

			lwork = -1;
			// eigvalue decomposition
			// find  the size of the memory required.
			dsyev_(GLOBAL::V, GLOBAL::U, &k, U, &k, D, &lworkopt, &lwork, &info);
			lwork = static_cast<integer> (lworkopt);
			work = new double[lwork];

			// compute eigvalue decomposition
			dsyev_(GLOBAL::V, GLOBAL::U, &k, U, &k, D, work, &lwork, &info);
			delete[] work;

			/*P_parallel = Q U */
			dgemm_(GLOBAL::N, GLOBAL::N, &n, &k, &k, &GLOBAL::DONE, Q, &n, U, &k, &GLOBAL::DZERO, P_parallel, &n);
		}

		if (k > 0)
		{
			/*Psig = R^T Q^T g, the g_parallel is used as a temporary data here*/
			dgemm_(GLOBAL::T, GLOBAL::N, &k, &GLOBAL::IONE, &n, &GLOBAL::DONE, Q, &n, const_cast<double *> (gf1->ObtainReadData()), &n, &GLOBAL::DZERO, g_parallel, &k);
			dgemm_(GLOBAL::N, GLOBAL::N, &k, &GLOBAL::IONE, &k, &GLOBAL::DONE, RT, &k, g_parallel, &k, &GLOBAL::DZERO, Psig, &k);
			/*g_parallel = P_parallel^T gf */
			dgemm_(GLOBAL::T, GLOBAL::N, &k, &GLOBAL::IONE, &n, &GLOBAL::DONE, P_parallel, &n, const_cast<double *> (gf1->ObtainReadData()), &n, &GLOBAL::DZERO, g_parallel, &k);
		}

		//ForDebug::Print("Psi:", Psi, n, k);//----
		//ForDebug::Print("PMGQ", PMGQ, k, k);//---
		//std::cout << "gamma:" << gamma << std::endl;//---
		//ForDebug::Print("gf1:", gf1->ObtainReadData(), gf1->Getlength());
		//std::cout << "delta:" << Delta << std::endl;//---

		for (integer i = 0; i < k; i++)
			Lambda[i] = D[i] + gamma;
		Lambda[k] = gamma;
		for (integer i = 0; i < k + 1; i++)
			if (fabs(Lambda[i]) <= tolLocalNewton)
				Lambda[i] = 0;
		double lambda_min = (Lambda[0] < gamma) ? Lambda[0] : gamma;
		double a_kp2 = std::sqrt(ddot_(&n, const_cast<double *> (gf1->ObtainReadData()), &GLOBAL::IONE, const_cast<double *> (gf1->ObtainReadData()), &GLOBAL::IONE) - ddot_(&k, g_parallel, &GLOBAL::IONE, g_parallel, &GLOBAL::IONE));
		if (a_kp2 * a_kp2 < tolLocalNewton)
			a_kp2 = 0;

		for (integer i = 0; i < k; i++)
			a_j[i] = g_parallel[i];
		a_j[k] = a_kp2;

		double tmpv = 0;
		for (integer i = 0; i < k + 1; i++)
			tmpv += a_j[i] * a_j[i] / Lambda[i] / Lambda[i];
		tmpv = sqrt(tmpv);

		double *pStar = eta2->ObtainWriteEntireData();

		double sigmaStar = 0;
		tCGstatus = L_TR_MIN;
		if (lambda_min > 0 && tmpv <= Delta)
		{
			ComputeSBySMW(gamma, Psig, RT); /*gf1, Psig, YMGS, PMGQ, RT*/
		}
		else
			if (lambda_min <= 0 && PhiBar_fg(-lambda_min, Delta, Lambda, a_j) > 0)
			{
				sigmaStar = -lambda_min;
				double *v = new double[k + 1];
				for (integer i = 0; i < k + 1; i++)
				{
					if (fabs(Lambda[i] + sigmaStar) > tolLocalNewton)
					{
						v[i] = a_j[i] / (Lambda[i] + sigmaStar);
					}
					else
					{
						v[i] = 0;
					}
				}
				dgemm_(GLOBAL::N, GLOBAL::N, &n, &GLOBAL::IONE, &k, &GLOBAL::DNONE, P_parallel, &n, v, &k, &GLOBAL::DZERO, pStar, &n);
				delete[] v;
				if (fabs(gamma + sigmaStar) > tolLocalNewton)
				{
					double *tmp = new double[n];
					dgemm_(GLOBAL::N, GLOBAL::N, &n, &GLOBAL::IONE, &k, &GLOBAL::DONE, P_parallel, &n, g_parallel, &k, &GLOBAL::DZERO, tmp, &n);
					double coef = 1.0 / (gamma + sigmaStar);
					daxpy_(&n, &coef, tmp, &GLOBAL::IONE, pStar, &GLOBAL::IONE);
					coef = -coef;
					daxpy_(&n, &coef, const_cast<double *> (gf1->ObtainReadData()), &GLOBAL::IONE, pStar, &GLOBAL::IONE);
					delete[] tmp;
				}

				if (lambda_min < 0)
				{
					double alpha = sqrt(Delta * Delta - ddot_(&n, pStar, &GLOBAL::IONE, pStar, &GLOBAL::IONE));

					double *zstar = new double[n];
					for (integer i = 0; i < n; i++)
						zstar[i] = 0;
					if (fabs(lambda_min - Lambda[0]) < tolLocalNewton)
					{
						double coef = alpha / sqrt(ddot_(&n, P_parallel, &GLOBAL::IONE, P_parallel, &GLOBAL::IONE));
						daxpy_(&n, &coef, P_parallel, &GLOBAL::IONE, zstar, &GLOBAL::IONE);
					}
					else
					{
						double norm_umin = 1;
						for (integer i = 0; i < k; i++)
						{
							dgemm_(GLOBAL::N, GLOBAL::T, &n, &GLOBAL::IONE, &k, &GLOBAL::DNONE, P_parallel, &n, P_parallel + i, &n, &GLOBAL::DZERO, zstar, &n);
							zstar[i] += 1;
							norm_umin = sqrt(ddot_(&n, zstar, &GLOBAL::IONE, zstar, &GLOBAL::IONE));
							if (norm_umin > tolLocalNewton)
								break;
						}
						double coef = alpha / norm_umin;
						dscal_(&n, &coef, zstar, &GLOBAL::IONE);
					}
					daxpy_(&n, &GLOBAL::DONE, zstar, &GLOBAL::IONE, pStar, &GLOBAL::IONE);
					delete[] zstar;
				}
				tCGstatus = TR_NEGCURVTURE;
			}
			else
			{
				if (lambda_min > 0)
				{
					sigmaStar = LocalNewton(0, 100, tolLocalNewton, Lambda, a_j);
				}
				else
				{
					double sigmaHat = 0;
					for (integer i = 0; i < k + 1; i++)
					{
						if (sigmaHat < a_j[i] / Delta - Lambda[i])
							sigmaHat = a_j[i] / Delta - Lambda[i];
					}
					if (sigmaHat > -lambda_min)
						sigmaStar = LocalNewton(sigmaHat, 100, tolLocalNewton, Lambda, a_j);
					else
						sigmaStar = LocalNewton(-lambda_min, 100, tolLocalNewton, Lambda, a_j);
				}
				ComputeSBySMW(gamma + sigmaStar, Psig, RT); /*gf1, Psig, YMGS, PMGQ, RT*/
				tCGstatus = TR_EXCREGION;
			}
	};

	void LRTRSR1::ComputeSBySMW(double tauStar, double *Psig, double * RT)
	{/*gamma, gf1, Psig, YMGS, PMGQ, RT*/
		integer n = gf1->Getlength(); /*number of rows*/
		integer k = Currentlength;    /*number of columns*/
		double *pStar = eta2->ObtainWriteEntireData();
		double coef = -1.0 / tauStar;

		if (k != 0)
		{
			double *vw = new double[k * k];
			integer *Ptmp = new integer[k];
			integer k2 = k * k;
			dgemm_(GLOBAL::N, GLOBAL::T, &k, &k, &k, &GLOBAL::DONE, RT, &k, RT, &k, &GLOBAL::DZERO, vw, &k);
			//ForDebug::Print("PsiPsi:", vw, k, k);//---
			dscal_(&k2, &tauStar, vw, &GLOBAL::IONE);
			double tauStar2 = tauStar * tauStar;
			daxpy_(&k2, &tauStar2, PMGQ, &GLOBAL::IONE, vw, &GLOBAL::IONE);
			//ForDebug::Print("vw:", vw, k, k);//---
			double *tmp = new double[k];
			integer info = 0;
			// LU decomposion for vw, vw = P * L * U, L and U are stored in vw, the permutation matrix is in Ptmp
			// details: http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html
			dgetrf_(&k, &k, vw, &k, Ptmp, &info);
			if (info < 0)
				printf("Error in dgetrf of LRTRSR1::ComputeSBySMW\n");

			dcopy_(&k, Psig, &GLOBAL::IONE, tmp, &GLOBAL::IONE);
			// solve linear system: vw * X = tmp using the LU decomposition results from dgetrf, then solution is stored in tmp.
			// details: http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f.html
			dgetrs_(GLOBAL::N, &k, &GLOBAL::IONE, vw, &k, Ptmp, tmp, &k, &info);
			if (info < 0)
				printf("Error in dgetrs of LRTRSR1::ComputeSBySMW\n");

			dgemm_(GLOBAL::N, GLOBAL::N, &n, &GLOBAL::IONE, &k, &GLOBAL::DONE, Psi, &n, tmp, &k, &GLOBAL::DZERO, pStar, &n);
			daxpy_(&n, &coef, const_cast<double *> (gf1->ObtainReadData()), &GLOBAL::IONE, pStar, &GLOBAL::IONE);

			delete[] vw;
			delete[] tmp;
			delete[] Ptmp;
			return;
		}
		gf1->CopyTo(eta2);
		dscal_(&n, &coef, pStar, &GLOBAL::IONE);
	};

	double LRTRSR1::PhiBar_fg(double nlambda_min, double Delta, double *Lambda, double *a_j, double *gf)
	{
		integer n = gf1->Getlength(); /*number of rows*/
		integer k = Currentlength;    /*number of columns*/
		double *tmp = new double[k + 1];
		for (integer i = 0; i < k + 1; i++)
			tmp[i] = Lambda[i] + nlambda_min;
		double eps_tol = 1e-10;
		double gradf = 0;

		bool flag = false;
		for (integer i = 0; i < k + 1; i++)
		{
			if (fabs(a_j[i]) < eps_tol || fabs(tmp[i]) < eps_tol)
			{
				flag = true;
				break;
			}
		}
		if (flag)
		{
			double pnorm2 = 0;
			for (integer i = 0; i < k + 1; i++)
			{
				if (fabs(a_j[i]) > eps_tol && fabs(tmp[i]) < eps_tol)
				{
					gradf = 1.0 / eps_tol;
					if (gf != nullptr)
						*gf = gradf;
					delete[] tmp;
					return -1.0 / Delta;
				}
				else
				if (fabs(a_j[i]) > eps_tol && fabs(tmp[i]) > eps_tol)
				{
					pnorm2 += a_j[i] * a_j[i] / tmp[i] / tmp[i];
					gradf += a_j[i] * a_j[i] / tmp[i] / tmp[i] / tmp[i];
				}
			}
			double norm = sqrt(pnorm2);
			gradf = gradf / norm / norm / norm;
			if (gf != nullptr)
				*gf = gradf;
			delete[] tmp;
			return 1.0 / norm - 1.0 / Delta;
		}

		for (integer i = 0; i < k + 1; i++)
			gradf += a_j[i] * a_j[i] / tmp[i] / tmp[i] / tmp[i];

		for (integer i = 0; i < k + 1; i++)
			tmp[i] = a_j[i] / tmp[i];
		integer length = k + 1;
		double norm = sqrt(ddot_(&length, tmp, &GLOBAL::IONE, tmp, &GLOBAL::IONE));

		gradf = gradf / norm / norm / norm;
		if (gf != nullptr)
			*gf = gradf;
		delete[] tmp;
		return 1.0 / norm - 1.0 / Delta;
	};

	double LRTRSR1::LocalNewton(double x0, integer maxIter, double tol, double *Lambda, double *a_j)
	{
		integer k = 0;
		double gf = 0;
		double fv = PhiBar_fg(x0, Delta, Lambda, a_j, &gf);
		while (fabs(fv) > std::numeric_limits<double>::epsilon() && k < maxIter)
		{
			x0 -= fv / gf;
			fv = PhiBar_fg(x0, Delta, Lambda, a_j, &gf);
			k++;
		}
		return x0;
	};

	void LRTRSR1::CheckParams(void)
	{
		SolversTR::CheckParams();
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		printf("LRTRSR1 METHOD PARAMETERS:\n");
		status = YES;
		printf("isconvex      :%15d[%s],\t", isconvex, status);
		status = (LengthSY >= 0) ? YES : NO;
		printf("LengthSY      :%15d[%s]\n", LengthSY, status);
	};

	void LRTRSR1::HessianEta(Vector *Eta, Vector *result)
	{
		HvLRTRSR1(Eta, result);
	};

	void LRTRSR1::UpdateData(void)
	{
		UpdateDataLRTRSR1();
	};

	void LRTRSR1::Acceptence(void)
	{
		for (integer i = 0; i < Currentlength; i++)
		{
			Mani->VectorTransport(x1, eta2, x2, S[i], S[i]);
			Mani->VectorTransport(x1, eta2, x2, Y[i], Y[i]);
		}
	};

	void LRTRSR1::PrintGenInfo(void)
	{
		Solvers::PrintGenInfo();
		printf("nH:%d,rho:%.2e,radius:%.3e,tCGstatus:%s,", nH, rho, Delta, tCGstatusSetnames[tCGstatus].c_str());
	};

	void LRTRSR1::PrintInfo(void)
	{
		printf("\n\tgamma:%.3e,inpss:%.3e,inpsy:%.3e,inpyy:%.3e,IsUpdateHessian:%d,", gamma, inpss, inpsy, inpyy, isupdated);
		printf("\n");
	};
}; /*end of ROPTLIB namespace*/
