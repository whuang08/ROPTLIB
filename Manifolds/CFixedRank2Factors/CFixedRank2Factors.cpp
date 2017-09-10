
#include "Manifolds/CFixedRank2Factors/CFixedRank2Factors.h"

/*Define the namespace*/
namespace ROPTLIB{

	CFixedRank2Factors::CFixedRank2Factors(integer inm, integer inn, integer inr) : ProductManifold(2,
		new Euclidean(2 * inm, inr), static_cast<integer> (1), new Euclidean(2 * inn, inr), static_cast<integer> (1))
	{
		m = inm;
		n = inn;
		r = inr;
		name.assign("Complex fixed-rank manifold by 2-factor representation");
		IsIntrApproach = true;
		delete EMPTYEXTR;
		delete EMPTYINTR;
		EMPTYEXTR = new CFR2Vector(m, n, r);
		EMPTYINTR = new CFR2Vector(m, n - r, r);
	};

	CFixedRank2Factors::~CFixedRank2Factors()
	{
		for (integer i = 0; i < numofmani; i++)
		{
			delete manifolds[i];
		}
	};

	void CFixedRank2Factors::ObtainIntr(Variable *x, Vector *etax, Vector *result) const
	{
		CFR2Variable *CFR2x = dynamic_cast<CFR2Variable *> (x);
		Variable* G = CFR2x->GetElement(0);
		Variable* H = CFR2x->GetElement(1);
		if (!G->TempDataExist("HHR"))
		{
			CpxNStQOrth::ComputeHHR(G);
		}
		if (!H->TempDataExist("HHR"))
		{
			CpxNStQOrth::ComputeHHR(H);
		}

		const double *etaxGTV = etax->ObtainReadData();
		const double *etaxHTV = etaxGTV + 2 * m * r;
		const SharedSpace *GHHR = G->ObtainReadTempData("HHR");
		const SharedSpace *GHHRTau = G->ObtainReadTempData("HHRTau");
		const SharedSpace *HHHR = H->ObtainReadTempData("HHR");
		const SharedSpace *HHHRTau = H->ObtainReadTempData("HHRTau");
		double *resultTV = result->ObtainWriteEntireData();
		const double *ptrGHHR = GHHR->ObtainReadData();
		const double *ptrGHHRTau = GHHRTau->ObtainReadData();
		const double *ptrHHHR = HHHR->ObtainReadData();
		const double *ptrHHHRTau = HHHRTau->ObtainReadData();

		/*for G component*/
		integer info;
		integer lwork = -1, Length = 2 * m * r;
		doublecomplex lworkopt;
		doublecomplex *Gtempspace = new doublecomplex[2 * m * r];
		// compute the size of space required in the dormqr
#ifndef MATLAB_MEX_FILE
		zunmqr_(GLOBAL::L, GLOBAL::C, &m, &r, &r, (doublecomplex *) (const_cast<double *> (ptrGHHR)), &m,
			(doublecomplex *)(const_cast<double *> (ptrGHHRTau)), Gtempspace, &m, &lworkopt, &lwork, &info);
#else
		zunmqr_(GLOBAL::L, GLOBAL::C, &m, &r, &r, (const_cast<double *> (ptrGHHR)), &m,
			(const_cast<double *> (ptrGHHRTau)), (double *)Gtempspace, &m, (double *)&lworkopt, &lwork, &info);
#endif
		lwork = static_cast<integer> (lworkopt.r);
		doublecomplex *work = new doublecomplex[lwork];
		// tempspace <- etaxTV, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		dcopy_(&Length, const_cast<double *> (etaxGTV), &GLOBAL::IONE, (double *)Gtempspace, &GLOBAL::IONE);
		// tempspace <- Q^T * tempspace, where Q is the orthogonal matrix defined as the product elementary reflectors defined by ptrHHR and ptrHHRTau,
		// details: http://www.netlib.org/lapack/explore-html/da/d82/dormqr_8f.html
#ifndef MATLAB_MEX_FILE
		zunmqr_(GLOBAL::L, GLOBAL::C, &m, &r, &r, (doublecomplex *)(const_cast<double *> (ptrGHHR)), &m,
			(doublecomplex *)(const_cast<double *> (ptrGHHRTau)), Gtempspace, &m, work, &lwork, &info);
#else
		zunmqr_(GLOBAL::L, GLOBAL::C, &m, &r, &r, (const_cast<double *> (ptrGHHR)), &m,
			(const_cast<double *> (ptrGHHRTau)), (double *)Gtempspace, &m, (double *)work, &lwork, &info);
#endif
		delete[] work;
		doublecomplex sign = { 0, 0 };
		for (integer i = 0; i < r; i++)
		{
			sign.r = (((doublecomplex *)ptrGHHR)[i + m * i].r >= 0) ? 1 : -1;
#ifndef MATLAB_MEX_FILE
			zscal_(&r, &sign, Gtempspace + i, &m);
#else
			zscal_(&r, (double *)&sign, (double *)(Gtempspace + i), &m);
#endif
		}

		doublecomplex *LG = new doublecomplex[r * r];
		for (integer i = 0; i < r; i++)
		{
			for (integer j = 0; j < i; j++)
			{
				LG[j + i * r].r = 0;
				LG[j + i * r].i = 0;
			}
			sign.r = (((doublecomplex *)ptrGHHR)[i + m * i].r >= 0) ? 1 : -1;
			for (integer j = i; j < r; j++)
			{
				LG[j + i * r].r = ((doublecomplex *)ptrGHHR)[i + m * j].r * sign.r;
				LG[j + i * r].i = ((doublecomplex *)ptrGHHR)[i + m * j].i * (-sign.r);
			}
		}

		/*for H component*/
		doublecomplex *Htempspace = new doublecomplex[2 * n * r];
		lwork = -1;
#ifndef MATLAB_MEX_FILE
		zunmqr_(GLOBAL::L, GLOBAL::C, &n, &r, &r, (doublecomplex *)(const_cast<double *> (ptrHHHR)), &n,
			(doublecomplex *)(const_cast<double *> (ptrHHHRTau)), Htempspace, &n, &lworkopt, &lwork, &info);
#else
		zunmqr_(GLOBAL::L, GLOBAL::C, &n, &r, &r, (const_cast<double *> (ptrHHHR)), &n,
			(const_cast<double *> (ptrHHHRTau)), (double *)Htempspace, &n, (double *)&lworkopt, &lwork, &info);
#endif
		lwork = static_cast<integer> (lworkopt.r);
		work = new doublecomplex[lwork];
		Length = 2 * n * r;
		// tempspace <- etaxTV, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		dcopy_(&Length, const_cast<double *> (etaxHTV), &GLOBAL::IONE, (double *)Htempspace, &GLOBAL::IONE);
		// tempspace <- Q^T * tempspace, where Q is the orthogonal matrix defined as the product elementary reflectors defined by ptrHHR and ptrHHRTau,
		// details: http://www.netlib.org/lapack/explore-html/da/d82/dormqr_8f.html
#ifndef MATLAB_MEX_FILE
		zunmqr_(GLOBAL::L, GLOBAL::C, &n, &r, &r, (doublecomplex *)(const_cast<double *> (ptrHHHR)), &n,
			(doublecomplex *)(const_cast<double *> (ptrHHHRTau)), Htempspace, &n, work, &lwork, &info);
#else
		zunmqr_(GLOBAL::L, GLOBAL::C, &n, &r, &r, (const_cast<double *> (ptrHHHR)), &n,
			(const_cast<double *> (ptrHHHRTau)), (double *)Htempspace, &n, (double *)work, &lwork, &info);
#endif
		delete[] work;
		sign.r = 0;
		sign.i = 0;
		for (integer i = 0; i < r; i++)
		{
			sign.r = (((doublecomplex *)ptrHHHR)[i + n * i].r >= 0) ? 1 : -1;
#ifndef MATLAB_MEX_FILE
			zscal_(&r, &sign, Htempspace + i, &n);
#else
			zscal_(&r, (double *)&sign, (double *)(Htempspace + i), &n);
#endif
		}

		doublecomplex *LH = new doublecomplex[r * r];
		for (integer i = 0; i < r; i++)
		{
			for (integer j = 0; j < i; j++)
			{
				LH[j + i * r].r = 0;
				LH[j + i * r].i = 0;
			}
			sign.r = (((doublecomplex *)ptrHHHR)[i + n * i].r >= 0) ? 1 : -1;
			for (integer j = i; j < r; j++)
			{
				LH[j + i * r].r = ((doublecomplex *)ptrHHHR)[i + n * j].r * sign.r;
				LH[j + i * r].i = ((doublecomplex *)ptrHHHR)[i + n * j].i * (-sign.r);
			}
		}

		/*Compute representations*/

		doublecomplex *GtempspaceL = new doublecomplex[m * r];

		Matrix GMtL((double *)GtempspaceL, m, r), MLH((double *)LH, r, r), MGtempspace((double *)Gtempspace, m, r);
		Matrix::CGEMM(GLOBAL::ZONE, MGtempspace, false, MLH, false, GLOBAL::ZZERO, GMtL);

		doublecomplex *HtempspaceL = new doublecomplex[n * r];

		Matrix HMtL((double *)HtempspaceL, n, r), MLG((double *)LG, r, r), MHtempspace((double *)Htempspace, n, r);
		Matrix::CGEMM(GLOBAL::ZONE, MHtempspace, false, MLG, false, GLOBAL::ZZERO, HMtL);
		delete[] LG;
		delete[] LH;
		delete[] Gtempspace;
		delete[] Htempspace;

		//ForDebug::Print("GtempspaceL", (double *)GtempspaceL, 2 * m, r);//---
		//ForDebug::Print("HtempspaceL", (double *)HtempspaceL, 2 * n, r);//---

		double r2 = sqrt(2.0);
		integer idx = 0;
		for (integer i = 0; i < r; i++)
		{
			for (integer j = 0; j < r; j++)
			{
				resultTV[idx] = (GtempspaceL[j + i * m].r + HtempspaceL[i + j * n].r) / r2;
				idx++;
				resultTV[idx] = (GtempspaceL[j + i * m].i - HtempspaceL[i + j * n].i) / r2;
				idx++;
			}
		}

		for (integer i = 0; i < r; i++)
		{
			for (integer j = r; j < m; j++)
			{
				resultTV[idx] = GtempspaceL[j + i * m].r;
				idx++;
				resultTV[idx] = GtempspaceL[j + i * m].i;
				idx++;
			}
		}
		for (integer i = 0; i < r; i++)
		{
			for (integer j = r; j < n; j++)
			{
				resultTV[idx] = HtempspaceL[j + i * n].r;
				idx++;
				resultTV[idx] = HtempspaceL[j + i * n].i;
				idx++;
			}
		}
		delete[] GtempspaceL;
		delete[] HtempspaceL;
	};

	void CFixedRank2Factors::ObtainExtr(Variable *x, Vector *intretax, Vector *result) const
	{
		CFR2Variable *CFR2x = dynamic_cast<CFR2Variable *> (x);
		Variable* G = CFR2x->GetElement(0);
		Variable* H = CFR2x->GetElement(1);
		if (!G->TempDataExist("HHR"))
		{
			CpxNStQOrth::ComputeHHR(G);
		}
		if (!H->TempDataExist("HHR"))
		{
			CpxNStQOrth::ComputeHHR(H);
		}

		const double *ptrintretax = intretax->ObtainReadData();
		const SharedSpace *GHHR = G->ObtainReadTempData("HHR");
		const SharedSpace *GHHRTau = G->ObtainReadTempData("HHRTau");
		const SharedSpace *HHHR = H->ObtainReadTempData("HHR");
		const SharedSpace *HHHRTau = H->ObtainReadTempData("HHRTau");
		doublecomplex *Gresult = (doublecomplex *) (result->ObtainWriteEntireData());
		doublecomplex *Hresult = Gresult + m * r;
		const double *ptrGHHR = GHHR->ObtainReadData();
		const double *ptrGHHRTau = GHHRTau->ObtainReadData();
		const double *ptrHHHR = HHHR->ObtainReadData();
		const double *ptrHHHRTau = HHHRTau->ObtainReadData();



		integer info;
		integer idx = 0;
		double r2 = sqrt(2.0);
		for (integer i = 0; i < r; i++)
		{
			for (integer j = 0; j < r; j++)
			{
				Gresult[j + i * m].r = ptrintretax[idx] / r2;
				Hresult[i + j * n].r = Gresult[j + i * m].r;
				idx++;
				Gresult[j + i * m].i = ptrintretax[idx] / r2;
				Hresult[i + j * n].i = -Gresult[j + i * m].i;
				idx++;
			}
		}

		for (integer i = 0; i < r; i++)
		{
			for (integer j = r; j < m; j++)
			{
				Gresult[j + i * m].r = ptrintretax[idx];
				idx++;
				Gresult[j + i * m].i = ptrintretax[idx];
				idx++;
			}
		}
		for (integer i = 0; i < r; i++)
		{
			for (integer j = r; j < n; j++)
			{
				Hresult[j + i * n].r = ptrintretax[idx];
				idx++;
				Hresult[j + i * n].i = ptrintretax[idx];
				idx++;
			}
		}

		//ForDebug::Print("Gresult", (double *)Gresult, 2 * m, r);//---
		//ForDebug::Print("Hresult", (double *)Hresult, 2 * n, r);//---

		/*for eta_G*/

		doublecomplex sign = { 0, 0 };
		for (integer i = 0; i < r; i++)
		{
			sign.r = (((doublecomplex *)ptrGHHR)[i + m * i].r >= 0) ? 1 : -1;
			// result(i, :) <- sign * result(i, :), details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
#ifndef MATLAB_MEX_FILE
			zscal_(&r, &sign, Gresult + i, &m);
#else
			zscal_(&r, (double *)&sign, ((double *) Gresult) + 2 * i, &m);
#endif
		}
		integer lwork = -1;
		doublecomplex lworkopt;
		// compute the size of space required in the dormqr
#ifndef MATLAB_MEX_FILE
		zunmqr_(GLOBAL::L, GLOBAL::N, &m, &r, &r, (doublecomplex *)(const_cast<double *> (ptrGHHR)), &m,
			(doublecomplex *)(const_cast<double *> (ptrGHHRTau)), Gresult, &m, &lworkopt, &lwork, &info);
#else
		zunmqr_(GLOBAL::L, GLOBAL::N, &m, &r, &r, (const_cast<double *> (ptrGHHR)), &m,
			(const_cast<double *> (ptrGHHRTau)), (double *) Gresult, &m, (double *)&lworkopt, &lwork, &info);
#endif
		lwork = static_cast<integer> (lworkopt.r);

		doublecomplex *work = new doublecomplex[lwork];
		// resultTV <- Q * resultTV, where Q is the orthogonal matrix defined as the product elementary reflectors defined by ptrHHR and ptrHHRTau,
		// details: http://www.netlib.org/lapack/explore-html/da/d82/dormqr_8f.html
#ifndef MATLAB_MEX_FILE
		zunmqr_(GLOBAL::L, GLOBAL::N, &m, &r, &r, (doublecomplex *)(const_cast<double *> (ptrGHHR)), &m,
			(doublecomplex *)(const_cast<double *> (ptrGHHRTau)), Gresult, &m, work, &lwork, &info);
#else
		zunmqr_(GLOBAL::L, GLOBAL::N, &m, &r, &r, (const_cast<double *> (ptrGHHR)), &m,
			(const_cast<double *> (ptrGHHRTau)), (double *) Gresult, &m, (double *)work, &lwork, &info);
#endif
		delete[] work;

		/*for eta_H*/
		sign.r = 0;
		sign.i = 0;
		for (integer i = 0; i < r; i++)
		{
			sign.r = (((doublecomplex *)ptrHHHR)[i + n * i].r >= 0) ? 1 : -1;
			// result(i, :) <- sign * result(i, :), details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
#ifndef MATLAB_MEX_FILE
			zscal_(&r, &sign, Hresult + i, &n);
#else
			zscal_(&r, (double *)&sign, ((double *)Hresult) + 2 * i, &n);
#endif
		}
		lwork = -1;
		// compute the size of space required in the dormqr
#ifndef MATLAB_MEX_FILE
		zunmqr_(GLOBAL::L, GLOBAL::N, &n, &r, &r, (doublecomplex *)(const_cast<double *> (ptrHHHR)), &n,
			(doublecomplex *)(const_cast<double *> (ptrHHHRTau)), Hresult, &n, &lworkopt, &lwork, &info);
#else
		zunmqr_(GLOBAL::L, GLOBAL::N, &n, &r, &r, (const_cast<double *> (ptrHHHR)), &n,
			(const_cast<double *> (ptrHHHRTau)), (double *)Hresult, &n, (double *)&lworkopt, &lwork, &info);
#endif
		lwork = static_cast<integer> (lworkopt.r);

		work = new doublecomplex[lwork];
		// resultTV <- Q * resultTV, where Q is the orthogonal matrix defined as the product elementary reflectors defined by ptrHHR and ptrHHRTau,
		// details: http://www.netlib.org/lapack/explore-html/da/d82/dormqr_8f.html
#ifndef MATLAB_MEX_FILE
		zunmqr_(GLOBAL::L, GLOBAL::N, &n, &r, &r, (doublecomplex *)(const_cast<double *> (ptrHHHR)), &n,
			(doublecomplex *)(const_cast<double *> (ptrHHHRTau)), Hresult, &n, work, &lwork, &info);
#else
		zunmqr_(GLOBAL::L, GLOBAL::N, &n, &r, &r, (const_cast<double *> (ptrHHHR)), &n,
			(const_cast<double *> (ptrHHHRTau)), (double *)Hresult, &n, (double *)work, &lwork, &info);
#endif
		delete[] work;

		doublecomplex *LG = new doublecomplex[2 * r * r + n * r + m * r];
		doublecomplex *LH = LG + r * r;
		doublecomplex *Gr_T = LH + r * r;
		doublecomplex *Hr_T = Gr_T + m * r;
		for (integer i = 0; i < r; i++)
		{
			for (integer j = 0; j < i; j++)
			{
				LG[j + i * r].r = 0;
				LG[j + i * r].i = 0;
			}
			sign.r = (((doublecomplex *)ptrGHHR)[i + m * i].r >= 0) ? 1 : -1;
			for (integer j = i; j < r; j++)
			{
				LG[j + i * r].r = ((doublecomplex *)ptrGHHR)[i + m * j].r * sign.r;
				LG[j + i * r].i = ((doublecomplex *)ptrGHHR)[i + m * j].i * (-sign.r);
			}
		}

		/*Gr_T <-  Gresult transpose conjugate*/
		for (integer i = 0; i < m; i++)
		{
			for (integer j = 0; j < r; j++)
			{
				Gr_T[j + i * r].r = ((doublecomplex *)Gresult)[i + j * m].r;
				Gr_T[j + i * r].i = -((doublecomplex *)Gresult)[i + j * m].i;
			}
		}

		for (integer i = 0; i < r; i++)
		{
			for (integer j = 0; j < i; j++)
			{
				LH[j + i * r].r = 0;
				LH[j + i * r].i = 0;
			}
			sign.r = (((doublecomplex *)ptrHHHR)[i + n * i].r >= 0) ? 1 : -1;
			for (integer j = i; j < r; j++)
			{
				LH[j + i * r].r = ((doublecomplex *)ptrHHHR)[i + n * j].r * sign.r;
				LH[j + i * r].i = ((doublecomplex *)ptrHHHR)[i + n * j].i * (-sign.r);
			}
		}

		/*Hr_T <-  Gresult transpose conjugate*/
		for (integer i = 0; i < n; i++)
		{
			for (integer j = 0; j < r; j++)
			{
				Hr_T[j + i * r].r = ((doublecomplex *)Hresult)[i + j * n].r;
				Hr_T[j + i * r].i = -((doublecomplex *)Hresult)[i + j * n].i;
			}
		}

#ifndef MATLAB_MEX_FILE
		/*solve linear system LH^H M = Gr_T, the solution M is stored in Gr_T*/
		ztrtrs_(GLOBAL::L, GLOBAL::C, GLOBAL::N, &r, &m, LH, &r, Gr_T, &r, &info);
#else
		ztrtrs_(GLOBAL::L, GLOBAL::C, GLOBAL::N, &r, &m, (double *)LH, &r, (double *)Gr_T, &r, &info);
#endif
		/*resultTV <-  r_T transpose conjugate*/
		for (integer i = 0; i < m; i++)
		{
			for (integer j = 0; j < r; j++)
			{
				((doublecomplex *)Gresult)[i + j * m].r = Gr_T[j + i * r].r;
				((doublecomplex *)Gresult)[i + j * m].i = -Gr_T[j + i * r].i;
			}
		}

#ifndef MATLAB_MEX_FILE
		/*solve linear system LH^H M = Gr_T, the solution M is stored in Gr_T*/
		ztrtrs_(GLOBAL::L, GLOBAL::C, GLOBAL::N, &r, &n, LG, &r, Hr_T, &r, &info);
#else
		ztrtrs_(GLOBAL::L, GLOBAL::C, GLOBAL::N, &r, &n, (double *)LG, &r, (double *)Hr_T, &r, &info);
#endif
		/*resultTV <-  r_T transpose conjugate*/
		for (integer i = 0; i < n; i++)
		{
			for (integer j = 0; j < r; j++)
			{
				((doublecomplex *)Hresult)[i + j * n].r = Hr_T[j + i * r].r;
				((doublecomplex *)Hresult)[i + j * n].i = -Hr_T[j + i * r].i;
			}
		}
		delete[] LG;

	};

	void CFixedRank2Factors::Retraction(Variable *x, Vector *etax, Variable *result, double stepsize) const
	{
		if (IsIntrApproach)
		{
			Vector *exetax = EMPTYEXTR->ConstructEmpty();
			ObtainExtr(x, etax, exetax);
			SetIsIntrApproach(false);
			ProductManifold::Retraction(x, exetax, result, stepsize);
			SetIsIntrApproach(true);
			delete exetax;
		}
		else
		{
			ProductManifold::Retraction(x, etax, result, stepsize);
		}
	};

	void CFixedRank2Factors::coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
	{
		xiy->CopyTo(result);
		printf("Warning: CFixedRank2Factors::coTangentVector has not been done!\n");
	};

	void CFixedRank2Factors::DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const
	{
		xix->CopyTo(result);
	};

	void CFixedRank2Factors::ExtrProjection(Variable *x, Vector *etax, Vector *result) const
	{
		CFR2Variable *CFR2x = dynamic_cast<CFR2Variable *> (x);
		Variable* G = CFR2x->GetElement(0);
		Variable* H = CFR2x->GetElement(1);
		const double *Gptr = G->ObtainReadData();
		const double *Hptr = H->ObtainReadData();
		const double *etaxGTV = etax->ObtainReadData();
		const double *etaxHTV = etaxGTV + 2 * m * r;
		

		/*Compute (H^*H)^{-1} H^* eta_H */
		doublecomplex *HH = new doublecomplex[4 * r * r];
		doublecomplex *HV = HH + r * r;
		doublecomplex *GG = HV + r * r;
		doublecomplex *GV = GG + r * r;
		Matrix MH(Hptr, n, r), MHH((double*)HH, r, r), MV(etaxHTV, n, r), MHV((double*)HV, r, r);
		// MHH <- MH^* MH
		Matrix::CGEMM(GLOBAL::ZONE, MH, true, MH, false, GLOBAL::ZZERO, MHH);
		// MHV <- MH^* MV
		Matrix::CGEMM(GLOBAL::ZONE, MH, true, MV, false, GLOBAL::ZZERO, MHV);

		integer info;
#ifndef MATLAB_MEX_FILE
		// solve for (MHH)^{-1} MHV
		zpotrf_(GLOBAL::L, &r, HH, &r, &info);
		zpotrs_(GLOBAL::L, &r, &r, HH, &r, HV, &r, &info);
#else
		// solve for (MHH)^{-1} MHV
		zpotrf_(GLOBAL::L, &r, (double *)HH, &r, &info);
		zpotrs_(GLOBAL::L, &r, &r, (double *)HH, &r, (double *)HV, &r, &info);
#endif
		if (info != 0)
		{
			printf("warning: zpotrs failed in CFixedRank2Factors::ExtrProjection with info:%d!\n", info);
		}

		/*Compute (G^*G)^{-1} G^* eta_G */
		Matrix MG(Gptr, m, r), MGG((double*)GG, r, r), MV2(etaxGTV, m, r), MGV((double*)GV, r, r);
		// MGG <- MG^* MG
		Matrix::CGEMM(GLOBAL::ZONE, MG, true, MG, false, GLOBAL::ZZERO, MGG);
		// MGV <- MG^* MV2
		Matrix::CGEMM(GLOBAL::ZONE, MG, true, MV2, false, GLOBAL::ZZERO, MGV);

#ifndef MATLAB_MEX_FILE
		// solve for (MGG)^{-1} MGV
		zpotrf_(GLOBAL::L, &r, GG, &r, &info);
		zpotrs_(GLOBAL::L, &r, &r, GG, &r, GV, &r, &info);
#else
		// solve for (MGG)^{-1} MGV
		zpotrf_(GLOBAL::L, &r, (double *)GG, &r, &info);
		zpotrs_(GLOBAL::L, &r, &r, (double *)GG, &r, (double *)GV, &r, &info);
#endif
		if (info != 0)
		{
			printf("warning: zpotrs failed in CFixedRank2Factors::ExtrProjection with info:%d!\n", info);
		}

		/*GV = (- HV^* + GV)/2 */
		for (integer i = 0; i < r; i++)
		{
			for (integer j = 0; j < r; j++)
			{
				GV[i + j * r].r -= HV[j + i * r].r;
				GV[i + j * r].r /= 2.0;
				GV[i + j * r].i += HV[j + i * r].i;
				GV[i + j * r].i /= 2.0;
			}
		}

		etax->CopyTo(result);
		double *resultTV = result->ObtainWritePartialData();
		Matrix MGresultTV(resultTV, m, r), MHresultTV(resultTV + 2 * m * r, n, r);
		Matrix::CGEMM(GLOBAL::ZNONE, MG, false, MGV, false, GLOBAL::ZONE, MGresultTV);
		Matrix::CGEMM(GLOBAL::ZONE, MH, false, MGV, true, GLOBAL::ZONE, MHresultTV);

		delete[] HH;
	};

	void CFixedRank2Factors::EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const
	{
		if (prob->GetUseHess())
		{
			Vector *segf = egf->ConstructEmpty();
			segf->NewMemoryOnWrite(); // I don't remember the reason. It seems to be required.
			egf->CopyTo(segf);
			SharedSpace *Sharedegf = new SharedSpace(segf);
			x->AddToTempData("EGrad", Sharedegf);
		}
		ExtrProjection(x, egf, gf);
	};

	void CFixedRank2Factors::EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const
	{
		ExtrProjection(x, exix, xix);
	};
}; /*end of ROPTLIB namespace*/
