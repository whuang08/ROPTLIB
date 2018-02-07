
#include "Problems/CFR2BlindDecon2D/CFR2BlindDecon2D.h"

#ifdef ROPTLIB_WITH_FFTW

/*Define the namespace*/
namespace ROPTLIB {

	CFR2BlindDecon2D::CFR2BlindDecon2D(double *iny, double *inB, integer innzmaxB, int *inirB, int *injcB, bool inisBsparse,
		double *inC, integer innzmaxC, int *inirC, int *injcC, bool inisCsparse, integer inn1, integer inn2, integer inr, double inrho, double ind, double inmu)
	{
		y = iny;

		B = inB;
		nzmaxB = innzmaxB;
		irB = inirB;
		jcB = injcB;
		isBsparse = inisBsparse;

		C = inC;
		nzmaxC = innzmaxC;
		irC = inirC;
		jcC = injcC;
		isCsparse = inisCsparse;

		n1 = inn1;
		n2 = inn2;
		L = n1 * n2;
		r = inr;

		rho = inrho;
		d = ind;
		mu = inmu;

		flags = FFTW_ESTIMATE;
		p = fftw_plan_dft_2d(n2, n1, nullptr, nullptr, FFTW_FORWARD, flags);

		if (isBsparse)
		{
			sB = BLAS_zuscr_begin(L, L);
			BLAS_zuscr_insert_entries(sB, nzmaxB, B, irB, jcB);
			BLAS_zuscr_end(sB);
		}
		if (isCsparse && C != nullptr)
		{
			sC = BLAS_zuscr_begin(L, L);
			BLAS_zuscr_insert_entries(sC, nzmaxC, C, irC, jcC);
			BLAS_zuscr_end(sC);
		}
		//log2L = static_cast<int> (log(static_cast<float> (L)) / log(static_cast<float> (2)) + 0.5);
		if (static_cast<int> (pow(2.0, static_cast<int>(log(static_cast<float> (L)) / log(static_cast<float> (2)) + 0.5))) - L != 0)
			printf("Warning: L must be a power of 2!\n");
	};

	CFR2BlindDecon2D::~CFR2BlindDecon2D(void)
	{
		fftw_destroy_plan(p);
		if (isBsparse)
			BLAS_usds(sB);
		if (isCsparse && C != nullptr)
			BLAS_usds(sC);
	};

	void CFR2BlindDecon2D::fftwrapper(integer inn1, integer inn2, fftw_complex *in, fftw_complex *out, int sign) const
	{
		if (r != 1)
			printf("Warning: CFR2BlindDecon2D::fftwrapper: suppose r = 1!\n");
		p = fftw_plan_dft_2d(inn2, inn1, in, out, sign, flags);
		fftw_execute(p);
	};

	void CFR2BlindDecon2D::BtimesU(double *U, double *result) const
	{
		if (isBsparse)
		{
			for (integer i = 0; i < r * 2 * L; i++)
				result[i] = 0;
			BLAS_zusmm(blas_colmajor, blas_no_trans, r, &GLOBAL::ZONE, sB, U, L, result, L);
		}
		else
		{
#ifndef MATLAB_MEX_FILE
			zgemm_(GLOBAL::N, GLOBAL::N, &L, &r, &L, &GLOBAL::ZONE, (doublecomplex*)B, &L, (doublecomplex*)U, &L, &GLOBAL::ZZERO, (doublecomplex *)result, &L);
#else
			zgemm_(GLOBAL::N, GLOBAL::N, &L, &r, &L, (double *)&GLOBAL::ZONE, B, &L, U, &L, (double *)&GLOBAL::ZZERO, result, &L);
#endif
		}
	};

	void CFR2BlindDecon2D::CtimesV(double *V, double *result) const
	{
		if (isCsparse)
		{
			for (integer i = 0; i < r * 2 * L; i++)
				result[i] = 0;
			BLAS_zusmm(blas_colmajor, blas_no_trans, r, &GLOBAL::ZONE, sC, V, L, result, L);
		}
		else
		{
#ifndef MATLAB_MEX_FILE
			zgemm_(GLOBAL::N, GLOBAL::N, &L, &r, &L, &GLOBAL::ZONE, (doublecomplex*)C, &L, (doublecomplex*)V, &L, &GLOBAL::ZZERO, (doublecomplex*)result, &L);
#else
			zgemm_(GLOBAL::N, GLOBAL::N, &L, &r, &L, (double *)&GLOBAL::ZONE, C, &L, V, &L, (double *)&GLOBAL::ZZERO, result, &L);
#endif
		}

		if (r != 1)
			printf("Warning: CFR2BlindDecon2D::CtimesV: suppose r = 1!\n");
		if (n1 * n2 != L)
			printf("Warning: n1 * n2 must equal L!\n");

		haarFWT_2d_inverse(n1, n2, (doublecomplex *)result);
		//haarFWT_1d_inverse(n1 * n2, (doublecomplex *)result);
	};

	void CFR2BlindDecon2D::BHtimesM(double *M, double *result) const
	{
		if (isBsparse)
		{
			for (integer i = 0; i < 2 * r * L; i++)
				result[i] = 0;
			BLAS_zusmm(blas_colmajor, blas_conj_trans, r, &GLOBAL::ZONE, sB, M, L, result, L);
		}
		else
		{
#ifndef MATLAB_MEX_FILE
			zgemm_(GLOBAL::C, GLOBAL::N, &L, &r, &L, &GLOBAL::ZONE, (doublecomplex *)B, &L, (doublecomplex *)M, &L, &GLOBAL::ZZERO, (doublecomplex *)result, &L);
#else
			zgemm_(GLOBAL::C, GLOBAL::N, &L, &r, &L, (double *)&GLOBAL::ZONE, B, &L, M, &L, (double *)&GLOBAL::ZZERO, result, &L);
#endif
		}
	};

	void CFR2BlindDecon2D::CHtimesM(double *M, double *result) const
	{ /*Data in M will be destroy in this function*/
		if (r != 1)
			printf("Warning: CFR2BlindDeconvolution::CHtimesM: suppose r = 1!\n");

		if (n1 * n2 != L)
			printf("Warning: n1 * n2 must equal L!\n");

		haarFWT_2d(n1, n2, (doublecomplex *)M);
		//haarFWT_1d(n1 * n2, (doublecomplex *)M);

		if (isCsparse)
		{
			for (integer i = 0; i < 2 * r * L; i++)
				result[i] = 0;
			BLAS_zusmm(blas_colmajor, blas_conj_trans, r, &GLOBAL::ZONE, sC, M, L, result, L);
		}
		else
		{
#ifndef MATLAB_MEX_FILE
			zgemm_(GLOBAL::C, GLOBAL::N, &L, &r, &L, &GLOBAL::ZONE, (doublecomplex *)C, &L, (doublecomplex *)M, &L, &GLOBAL::ZZERO, (doublecomplex *)result, &L);
#else
			zgemm_(GLOBAL::C, GLOBAL::N, &L, &r, &L, (double *)&GLOBAL::ZONE, C, &L, M, &L, (double *)&GLOBAL::ZZERO, result, &L);
#endif
		}
	};

	double CFR2BlindDecon2D::f(Variable *x) const
	{
		const double *xptr = x->ObtainReadData();
		const double *Uptr = xptr;
		const double *Vptr = xptr + 2 * L * r;
		SharedSpace *BUCV = new SharedSpace(1, L * r * 4);
		double *BU = BUCV->ObtainWriteEntireData();
		double *CV = BU + L * r * 2;
		integer LL2 = 2 * L;
		BtimesU(const_cast<double *>(Uptr), BU);
		fftwrapper(n1, n2, (fftw_complex *)BU, (fftw_complex *)BU, FFTW_FORWARD);
		CtimesV(const_cast<double *> (Vptr), CV);
		fftwrapper(n1, n2, (fftw_complex *)CV, (fftw_complex *)CV, FFTW_BACKWARD);

		SharedSpace *ymdiagBXA = new SharedSpace(1, L * 2);
		double *tmpptr = ymdiagBXA->ObtainWriteEntireData();
		for (integer i = 0; i < L; i++) // this for loop can be improved
		{
			tmpptr[2 * i] = y[2 * i] - (BU[2 * i] * CV[2 * i] + BU[2 * i + 1] * CV[2 * i + 1]);
			tmpptr[2 * i + 1] = y[2 * i + 1] - (BU[2 * i + 1] * CV[2 * i] - BU[2 * i] * CV[2 * i + 1]);
		}

		double penalty = 0;
		if (rho != 0)
		{
			double tmp = static_cast<double> (L) / 8 / d / d / mu / mu;
			double tmp2 = 0;
			SharedSpace *rownorm2BU = new SharedSpace(1, L + 1);
			double *rownorm2BUptr = rownorm2BU->ObtainWriteEntireData();
			integer length = r;
			for (integer i = 0; i < L; i++)
			{
				rownorm2BUptr[i] = ddot_(&length, BU + i * 2, &LL2, BU + i * 2, &LL2);
				rownorm2BUptr[i] += ddot_(&length, BU + i * 2 + 1, &LL2, BU + i * 2 + 1, &LL2);

			}
			length = L * 2 * r;
			rownorm2BUptr[L] = ddot_(&length, const_cast<double *> (Vptr), &GLOBAL::IONE, const_cast<double *> (Vptr), &GLOBAL::IONE);
			for (integer i = 0; i < L; i++)
			{
				tmp2 = ((tmp * rownorm2BUptr[i] * rownorm2BUptr[L] - 1 < 0) ? 0 : (tmp * rownorm2BUptr[i] * rownorm2BUptr[L] - 1));
				penalty += rho * tmp2 * tmp2;
			}
			x->AddToTempData("rownorm2BU", rownorm2BU);
		}
		double result = ddot_(&LL2, tmpptr, &GLOBAL::IONE, tmpptr, &GLOBAL::IONE) + penalty;

		if (UseGrad)
		{
			x->AddToTempData("ymdiagBXA", ymdiagBXA);
			x->AddToTempData("BUCV", BUCV);
		}
		else
		{
			delete ymdiagBXA;
			delete BUCV;
		}
		return result;
	};

	void CFR2BlindDecon2D::EucGrad(Variable *x, Vector *gf) const
	{
		const double *xptr = x->ObtainReadData();
		const double *Uptr = xptr;
		const double *Vptr = xptr + 2 * L * r;

		const SharedSpace *ymdiagBXA = x->ObtainReadTempData("ymdiagBXA");
		const double *tmpptr = ymdiagBXA->ObtainReadData();
		const SharedSpace *BUCV = x->ObtainReadTempData("BUCV");
		const double *BU = BUCV->ObtainReadData();
		const double *CV = BU + 2 * r * L;

		integer LL2 = L * 2;

		double *EGFV = new double[2 * L * r + 2 * L * r];
		double *EGFTU = EGFV + 2 * L * r;

		double *tmpp = new double[2 * r * L];
		integer length = 2 * r * L;

		/*Computing EFGV*/
		dcopy_(&length, const_cast<double *> (CV), &GLOBAL::IONE, tmpp, &GLOBAL::IONE);

#ifndef MATLAB_MEX_FILE
		for (integer i = 0; i < L; i++)
			zscal_(&r, (doublecomplex *)(tmpptr + 2 * i), (doublecomplex *)(tmpp + 2 * i), &L);
#else
		for (integer i = 0; i < L; i++)
			zscal_(&r, const_cast<double *> (tmpptr + 2 * i), tmpp + 2 * i, &L);
#endif
		fftwrapper(n1, n2, (fftw_complex *)tmpp, (fftw_complex *)tmpp, FFTW_BACKWARD);
		BHtimesM(tmpp, EGFV);

		/*Computing EFGTU*/
		dcopy_(&length, const_cast<double *> (BU), &GLOBAL::IONE, tmpp, &GLOBAL::IONE);
		dscal_(&L, &GLOBAL::DNONE, const_cast<double *> (tmpptr + 1), &GLOBAL::ITWO);
#ifndef MATLAB_MEX_FILE
		for (integer i = 0; i < L; i++)
			zscal_(&r, (doublecomplex *)(tmpptr + 2 * i), (doublecomplex *)(tmpp + 2 * i), &L);
#else
		for (integer i = 0; i < L; i++)
			zscal_(&r, const_cast<double *> (tmpptr + 2 * i), tmpp + 2 * i, &L);
#endif
		dscal_(&L, &GLOBAL::DNONE, const_cast<double *> (tmpptr + 1), &GLOBAL::ITWO);
		fftwrapper(n1, n2, (fftw_complex *)tmpp, (fftw_complex *)tmpp, FFTW_FORWARD);
		CHtimesM(tmpp, EGFTU);

		delete[] tmpp;

		length = 2 * r * L;
		double ntwo = -2;
		dscal_(&length, &ntwo, EGFV, &GLOBAL::IONE); /*gfh*/
		length = 2 * r * L;
		dscal_(&length, &ntwo, EGFTU, &GLOBAL::IONE); /*gfm*/

		if (rho != 0)
		{
			const SharedSpace *rownorm2BU = x->ObtainReadTempData("rownorm2BU");
			const double *rownorm2BUptr = rownorm2BU->ObtainReadData();

			double *gGh = new double[2 * L * r + 2 * L * r];
			double *DFBU = gGh + 2 * L * r;
			length = 2 * L * r;
			dcopy_(&length, const_cast<double *> (BU), &GLOBAL::IONE, DFBU, &GLOBAL::IONE);
			doublecomplex coef = { 0, 0 };
			double tmp = static_cast<double> (L) / 8 / d / d / mu / mu;

			for (integer i = 0; i < L; i++)
			{
				coef.r = 2.0 * tmp * rho * 2 * ((tmp * rownorm2BUptr[i] * rownorm2BUptr[L] - 1 < 0) ? 0 : tmp * rownorm2BUptr[i] * rownorm2BUptr[L] - 1) * rownorm2BUptr[L];
#ifndef MATLAB_MEX_FILE
				zscal_(&r, &coef, (doublecomplex *)(DFBU + 2 * r * i), &L);
#else
				zscal_(&r, (double *)(&coef), DFBU + 2 * r * i, &L);
#endif
			}
			fftwrapper(n1, n2, (fftw_complex *)DFBU, (fftw_complex *)DFBU, FFTW_BACKWARD);
			BHtimesM(DFBU, gGh);
			length = 2 * L * r;
			daxpy_(&length, &GLOBAL::DONE, gGh, &GLOBAL::IONE, EGFV, &GLOBAL::IONE);
			delete[] gGh;

			double coef2 = 0;
			for (integer i = 0; i < L; i++)
			{
				coef2 += 2.0 * tmp * rho * rownorm2BUptr[i] * 2 * ((tmp * rownorm2BUptr[i] * rownorm2BUptr[L] - 1 < 0) ? 0 : tmp * rownorm2BUptr[i] * rownorm2BUptr[L] - 1);
			}
			length = 2 * r * L;
			daxpy_(&length, &coef2, const_cast<double *> (Vptr), &GLOBAL::IONE, EGFTU, &GLOBAL::IONE);
		}

		/*Compute (H^*H)^{-1} H^* eta_H */
		SharedSpace *invVVUU = new SharedSpace(1, 4 * r * r);
		doublecomplex *VV = (doublecomplex *)(invVVUU->ObtainWriteEntireData());
		doublecomplex *UU = VV + r * r;

		Matrix MVV((double*)VV, r, r), MV(Vptr, L, r);
		// MVV <- MV^* MV
		Matrix::CGEMM(GLOBAL::ZONE, MV, true, MV, false, GLOBAL::ZZERO, MVV);
		Matrix MUU((double*)UU, r, r), MU(Uptr, L, r);
		// MUU <- MU^* MU
		Matrix::CGEMM(GLOBAL::ZONE, MU, true, MU, false, GLOBAL::ZZERO, MUU);

		double *gfh = gf->ObtainWriteEntireData();
		double *gfm = gfh + 2 * L * r;

		integer info = 0;
#ifndef MATLAB_MEX_FILE
		zpotrf_(GLOBAL::L, &r, VV, &r, &info);
		zpotri_(GLOBAL::L, &r, VV, &r, &info);
#else
		zpotrf_(GLOBAL::L, &r, (double *)VV, &r, &info);
		zpotri_(GLOBAL::L, &r, (double *)VV, &r, &info);
		//zpotrs_(GLOBAL::L, &rr, &KK, (double *)VV, &rr, (double *)EGFV, &rr, &info);
#endif
		Matrix MEGFV(EGFV, L, r), Mgfh(gfh, L, r);
		Matrix::CGEMM(GLOBAL::ZONE, MEGFV, false, MVV, false, GLOBAL::ZZERO, Mgfh);
		if (info != 0)
		{
			printf("warning: zpotrs failed in CFR2BlindDeconvolution::EucGrad with info:%d!\n", info);
		}

#ifndef MATLAB_MEX_FILE
		// solve for EGFTU (MUU)^{-1}
		zpotrf_(GLOBAL::L, &r, UU, &r, &info);
		zpotri_(GLOBAL::L, &r, UU, &r, &info);
#else
		// solve for EGFTU (MUU)^{-1}
		zpotrf_(GLOBAL::L, &r, (double *)UU, &r, &info);
		zpotri_(GLOBAL::L, &r, (double *)UU, &r, &info);
#endif
		Matrix MEGFTU(EGFTU, L, r), Mgfm(gfm, L, r);
		Matrix::CGEMM(GLOBAL::ZONE, MEGFTU, false, MUU, false, GLOBAL::ZZERO, Mgfm);
		if (info != 0)
		{
			printf("warning: zpotrs failed in CFR2BlindDeconvolution::EucGrad with info:%d!\n", info);
		}
		delete[] EGFV;
		if (UseHess)
		{
			x->AddToTempData("invVVUU", invVVUU);
		}
		else
		{
			delete invVVUU;
		}
	};

	void CFR2BlindDecon2D::EucHessianEta(Variable *x, Vector *etax, Vector *xix) const
	{ /*The Hessian of the penalty term is not considered here*/
		const double *xptr = x->ObtainReadData();
		const double *Uptr = xptr;
		const double *Vptr = xptr + 2 * L * r;
		const double *dUptr = etax->ObtainReadData();
		const double *dVptr = dUptr + 2 * L * r;

		const SharedSpace *BUCV = x->ObtainReadTempData("BUCV");
		const double *BU = BUCV->ObtainReadData();
		const double *CV = BU + 2 * r * L;

		double *BdU = new double[L * r * 4];
		double *CdV = BdU + L * r * 2;
		double *mdiagBetaxA = new double[L * 2];

		for (integer i = 0; i < 2 * L; i++)
		{
			mdiagBetaxA[i] = 0;
		}

		integer LL2 = 2 * L;
		BtimesU(const_cast<double *> (dUptr), BdU);
		fftwrapper(n1, n2, (fftw_complex *)BdU, (fftw_complex *)BdU, FFTW_FORWARD);
		for (integer i = 0; i < L; i++) // this for loop can be improved
		{
			mdiagBetaxA[2 * i] -= (BdU[2 * i] * CV[2 * i] + BdU[2 * i + 1] * CV[2 * i + 1]);
			mdiagBetaxA[2 * i + 1] -= (BdU[2 * i + 1] * CV[2 * i] - BdU[2 * i] * CV[2 * i + 1]);
		}

		CtimesV(const_cast<double *> (dVptr), CdV);
		fftwrapper(n1, n2, (fftw_complex *)CdV, (fftw_complex *)CdV, FFTW_BACKWARD);
		for (integer i = 0; i < L; i++) // this for loop can be improved
		{
			mdiagBetaxA[2 * i] -= (BU[2 * i] * CdV[2 * i] + BU[2 * i + 1] * CdV[2 * i + 1]);
			mdiagBetaxA[2 * i + 1] -= (BU[2 * i + 1] * CdV[2 * i] - BU[2 * i] * CdV[2 * i + 1]);
		}

		double *EHFV = new double[2 * L * r + 2 * L * r], *EHFTU = EHFV + 2 * L * r;
		double *tmpp = new double[2 * r * L];
		integer length = 2 * r * L;
		/*Computing EHFV*/
		dcopy_(&length, const_cast<double *> (CV), &GLOBAL::IONE, tmpp, &GLOBAL::IONE);

#ifndef MATLAB_MEX_FILE
		for (integer i = 0; i < L; i++)
			zscal_(&r, (doublecomplex *)(mdiagBetaxA + 2 * i), (doublecomplex *)(tmpp + 2 * i), &L);
#else
		for (integer i = 0; i < L; i++)
			zscal_(&r, mdiagBetaxA + 2 * i, tmpp + 2 * i, &L);
#endif
		fftwrapper(n1, n2, (fftw_complex *)tmpp, (fftw_complex *)tmpp, FFTW_BACKWARD);
		BHtimesM(tmpp, EHFV);

		/*Computing EHFTU*/
		dcopy_(&length, const_cast<double *> (BU), &GLOBAL::IONE, tmpp, &GLOBAL::IONE);
		dscal_(&L, &GLOBAL::DNONE, const_cast<double *> (mdiagBetaxA + 1), &GLOBAL::ITWO);
#ifndef MATLAB_MEX_FILE
		for (integer i = 0; i < L; i++)
			zscal_(&r, (doublecomplex *)(mdiagBetaxA + 2 * i), (doublecomplex *)(tmpp + 2 * i), &L);
#else
		for (integer i = 0; i < L; i++)
			zscal_(&r, mdiagBetaxA + 2 * i, tmpp + 2 * i, &L);
#endif
		dscal_(&L, &GLOBAL::DNONE, const_cast<double *> (mdiagBetaxA + 1), &GLOBAL::ITWO);
		fftwrapper(n1, n2, (fftw_complex *)tmpp, (fftw_complex *)tmpp, FFTW_FORWARD);
		CHtimesM(tmpp, EHFTU);

		delete[] tmpp;
		delete[] mdiagBetaxA;

		length = 2 * r * L;
		double ntwo = -2;
		dscal_(&length, &ntwo, EHFV, &GLOBAL::IONE);
		length = 2 * r * L;
		dscal_(&length, &ntwo, EHFTU, &GLOBAL::IONE);

		const SharedSpace *ymdiagBXA = x->ObtainReadTempData("ymdiagBXA");
		const double *tmpptr = ymdiagBXA->ObtainReadData();
		double *EGFdV = new double[2 * L * r], *EGFTdU = new double[2 * L * r];
		tmpp = new double[2 * r * L];
		length = 2 * r * L;
		/*Computing EGFdV*/
		dcopy_(&length, const_cast<double *> (CdV), &GLOBAL::IONE, tmpp, &GLOBAL::IONE);

#ifndef MATLAB_MEX_FILE
		for (integer i = 0; i < L; i++)
			zscal_(&r, (doublecomplex *)(tmpptr + 2 * i), (doublecomplex *)(tmpp + 2 * i), &L);
#else
		for (integer i = 0; i < L; i++)
			zscal_(&r, const_cast<double *> (tmpptr + 2 * i), tmpp + 2 * i, &L);
#endif
		fftwrapper(n1, n2, (fftw_complex *)tmpp, (fftw_complex *)tmpp, FFTW_BACKWARD);
		BHtimesM(tmpp, EGFdV);

		/*Computing EGFTdU*/
		dcopy_(&length, const_cast<double *> (BdU), &GLOBAL::IONE, tmpp, &GLOBAL::IONE);
		dscal_(&L, &GLOBAL::DNONE, const_cast<double *> (tmpptr + 1), &GLOBAL::ITWO);
#ifndef MATLAB_MEX_FILE
		for (integer i = 0; i < L; i++)
			zscal_(&r, (doublecomplex *)(tmpptr + 2 * i), (doublecomplex *)(tmpp + 2 * i), &L);
#else
		for (integer i = 0; i < L; i++)
			zscal_(&r, const_cast<double *> (tmpptr + 2 * i), tmpp + 2 * i, &L);
#endif
		dscal_(&L, &GLOBAL::DNONE, const_cast<double *> (tmpptr + 1), &GLOBAL::ITWO);
		fftwrapper(n1, n2, (fftw_complex *)tmpp, (fftw_complex *)tmpp, FFTW_FORWARD);
		CHtimesM(tmpp, EGFTdU);

		delete[] tmpp;

		ntwo = -2;
		length = 2 * r * L;
		daxpy_(&length, &ntwo, EGFdV, &GLOBAL::IONE, EHFV, &GLOBAL::IONE);
		length = 2 * r * L;
		daxpy_(&length, &ntwo, EGFTdU, &GLOBAL::IONE, EHFTU, &GLOBAL::IONE);

		delete[] EGFdV;
		delete[] EGFTdU;
		delete[] BdU;

		/*Compute (H^*H)^{-1} H^* eta_H */
		const SharedSpace *invVVUU = x->ObtainReadTempData("invVVUU");
		const doublecomplex *VV = (doublecomplex *)(invVVUU->ObtainReadData());
		const doublecomplex *UU = VV + r * r;

		Matrix MVV((double*)VV, r, r), MUU((double*)UU, r, r);
		double *Hfh = xix->ObtainWriteEntireData();
		double *Hfm = Hfh + 2 * L * r;

		Matrix MEHFV(EHFV, L, r), MHfh(Hfh, L, r);
		Matrix::CGEMM(GLOBAL::ZONE, MEHFV, false, MVV, false, GLOBAL::ZZERO, MHfh);

		Matrix MEHFTU(EHFTU, L, r), MHfm(Hfm, L, r);
		Matrix::CGEMM(GLOBAL::ZONE, MEHFTU, false, MUU, false, GLOBAL::ZZERO, MHfm);

		delete[] EHFV;
	};

}; /*end of ROPTLIB namespace*/
#endif
