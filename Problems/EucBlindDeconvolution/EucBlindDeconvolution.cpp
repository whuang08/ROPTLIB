
#include "Problems/EucBlindDeconvolution/EucBlindDeconvolution.h"

#ifdef ROPTLIB_WITH_FFTW

/*Define the namespace*/
namespace ROPTLIB {

	EucBlindDeconvolution::EucBlindDeconvolution(double *iny, double *inB, integer innzmaxB, int *inirB, int *injcB, bool inisBsparse,
		double *inC, integer innzmaxC, int *inirC, int *injcC, bool inisCsparse, integer inL, integer inK, integer inN, integer inr)
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

		L = inL;
		K = inK;
		N = inN;
		r = inr;
		flags = FFTW_ESTIMATE;
		p = fftw_plan_dft_1d(L, nullptr, nullptr, FFTW_FORWARD, flags);

		if (isBsparse)
		{
			sB = BLAS_duscr_begin(2 * L, K);
			BLAS_duscr_insert_entries(sB, nzmaxB, B, irB, jcB);
			BLAS_duscr_end(sB);
		}
		if (isCsparse)
		{
			sC = BLAS_duscr_begin(2 * L, N);
			BLAS_duscr_insert_entries(sC, nzmaxC, C, irC, jcC);
			BLAS_duscr_end(sC);
		}
		//double *testarr = new double[1024 * 1024 * 2];
		//for (integer i = 0; i < 1024 * 1024 * 2; i++)
		//	testarr[i] = 1;
		//unsigned long starttime = getTickCount();
		//fftwrapper(1024 * 1024, 1, (fftw_complex *)testarr, (fftw_complex *)testarr, FFTW_FORWARD);
		//printf("computational time:%e\n", static_cast<double>(getTickCount() - starttime) / CLK_PS);//---
		//std::cout << "=========================" << std::endl;//---
		//delete[] testarr;
	};

	EucBlindDeconvolution::~EucBlindDeconvolution(void)
	{
		fftw_destroy_plan(p);
		if (isBsparse)
			BLAS_usds(sB);
		if (isCsparse)
			BLAS_usds(sC);
	};

	void EucBlindDeconvolution::fftwrapper(integer inn, integer inr, fftw_complex *in, fftw_complex *out, int sign) const
	{
		for (integer i = 0; i < inr; i++)
		{
			p = fftw_plan_dft_1d(inn, in + inn * i, out + inn * i, sign, flags);
			fftw_execute(p);
		}
	};

	void EucBlindDeconvolution::BtimesU(double *U, double *result) const
	{
		if (isBsparse)
		{
			for (integer i = 0; i < r * 2 * L; i++)
				result[i] = 0;
			BLAS_dusmm(blas_colmajor, blas_no_trans, r, 1, sB, U, K, result, 2 * L);
		}
		else
		{
			integer LL2 = L * 2, rr = r, KK = K;
			dgemm_(GLOBAL::N, GLOBAL::N, &LL2, &rr, &KK, &GLOBAL::DONE, B, &LL2, U, &KK, &GLOBAL::DZERO, result, &LL2);
		}
	};

	void EucBlindDeconvolution::CtimesV(double *V, double *result) const
	{
		if (isCsparse)
		{
			for (integer i = 0; i < r * 2 * L; i++)
				result[i] = 0;
			BLAS_dusmm(blas_colmajor, blas_no_trans, r, 1, sC, V, N, result, 2 * L);
		}
		else
		{
			integer LL2 = L * 2, rr = r, NN = N;
			dgemm_(GLOBAL::N, GLOBAL::N, &LL2, &rr, &NN, &GLOBAL::DONE, C, &LL2, V, &NN, &GLOBAL::DZERO, result, &LL2);
		}
	};

	void EucBlindDeconvolution::BTtimesM(double *M, double *result) const
	{
		if (isBsparse)
		{
			for (integer i = 0; i < r * K; i++)
				result[i] = 0;
			BLAS_dusmm(blas_colmajor, blas_trans, r, 1, sB, M, 2 * L, result, K);
		}
		else
		{
			integer LL2 = L * 2, rr = r, KK = K;
			dgemm_(GLOBAL::T, GLOBAL::N, &KK, &rr, &LL2, &GLOBAL::DONE, B, &LL2, M, &LL2, &GLOBAL::DZERO, result, &KK);
		}
	};

	void EucBlindDeconvolution::CTtimesM(double *M, double *result) const
	{
		if (isCsparse)
		{
			for (integer i = 0; i < r * N; i++)
				result[i] = 0;
			BLAS_dusmm(blas_colmajor, blas_trans, r, 1, sC, M, 2 * L, result, N);
		}
		else
		{
			integer LL2 = L * 2, rr = r, NN = N;
			dgemm_(GLOBAL::T, GLOBAL::N, &NN, &rr, &LL2, &GLOBAL::DONE, C, &LL2, M, &LL2, &GLOBAL::DZERO, result, &NN);
		}
	};

	double EucBlindDeconvolution::f(Variable *x) const
	{/*x is a K by N matrix*/
		const double *xptr = x->ObtainReadData();
		const double *Uptr = xptr;
		const double *Vptr = xptr + K * r;
		double *BU = new double[L * r * 4];
		double *CV = BU + L * r * 2;
		integer LL2 = 2 * L, KK = K, NN = N, rr = r;
		BtimesU(const_cast<double *>(Uptr), BU);
		fftwrapper(L, r, (fftw_complex *)BU, (fftw_complex *)BU, FFTW_FORWARD);
		CtimesV(const_cast<double *> (Vptr), CV);
		fftwrapper(L, r, (fftw_complex *)CV, (fftw_complex *)CV, FFTW_FORWARD);

		SharedSpace *ymdiagBXA = new SharedSpace(1, L * 2);
		double *tmpptr = ymdiagBXA->ObtainWriteEntireData();
		//double *tmp = new double[4];
		for (integer i = 0; i < L; i++) // this for loop can be improved
		{
			tmpptr[2 * i] = y[2 * i] - (BU[2 * i] * CV[2 * i] - BU[2 * i + 1] * CV[2 * i + 1]);
			tmpptr[2 * i + 1] = y[2 * i + 1] - (BU[2 * i + 1] * CV[2 * i] + BU[2 * i] * CV[2 * i + 1]);
		}
		delete[] BU;
		double result = ddot_(&LL2, tmpptr, &GLOBAL::IONE, tmpptr, &GLOBAL::IONE);

		if (UseGrad)
		{
			x->AddToTempData("ymdiagBXA", ymdiagBXA);
		}
		else
		{
			delete ymdiagBXA;
		}
		return result;
	};

	void EucBlindDeconvolution::EucGrad(Variable *x, Vector *gf) const
	{
		const double *xptr = x->ObtainReadData();
		const double *Uptr = xptr;
		const double *Vptr = xptr + K * r;

		const SharedSpace *ymdiagBXA = x->ObtainReadTempData("ymdiagBXA");
		const double *tmpptr = ymdiagBXA->ObtainReadData();

		integer NN = N, LL = L, rr = r, KK = K, LL2 = L * 2;

		double *EGFV = gf->ObtainWriteEntireData();
		double *EGFTU = EGFV + K * r;

		double *tmpp = new double[2 * r * L];

		/*Computing EFGV*/

		CtimesV(const_cast<double *> (Vptr), tmpp);
		fftwrapper(L, r, (fftw_complex *)tmpp, (fftw_complex *)tmpp, FFTW_FORWARD);
		integer length = L * r;
		dscal_(&length, &GLOBAL::DNONE, tmpp + 1, &GLOBAL::ITWO);
#ifndef MATLAB_MEX_FILE
		for (integer i = 0; i < L; i++)
			zscal_(&rr, (doublecomplex *)(tmpptr + 2 * i), (doublecomplex *)(tmpp + 2 * i), &LL);
#else
		for (integer i = 0; i < L; i++)
			zscal_(&rr, const_cast<double *> (tmpptr + 2 * i), tmpp + 2 * i, &LL);
#endif
		fftwrapper(L, r, (fftw_complex *)tmpp, (fftw_complex *)tmpp, FFTW_BACKWARD);
		BTtimesM(tmpp, EGFV);

		/*Computing EFGTU*/
		BtimesU(const_cast<double *> (Uptr), tmpp);
		fftwrapper(L, r, (fftw_complex *)tmpp, (fftw_complex *)tmpp, FFTW_FORWARD);
		length = L * r;
		dscal_(&length, &GLOBAL::DNONE, tmpp + 1, &GLOBAL::ITWO);
#ifndef MATLAB_MEX_FILE
		for (integer i = 0; i < L; i++)
			zscal_(&rr, (doublecomplex *)(tmpptr + 2 * i), (doublecomplex *)(tmpp + 2 * i), &LL);
#else
		for (integer i = 0; i < L; i++)
			zscal_(&rr, const_cast<double *> (tmpptr + 2 * i), tmpp + 2 * i, &LL);
#endif
		fftwrapper(L, r, (fftw_complex *)tmpp, (fftw_complex *)tmpp, FFTW_BACKWARD);
		CTtimesM(tmpp, EGFTU);

		delete[] tmpp;

		length = r * KK;
		double ntwo = -2;
		dscal_(&length, &ntwo, EGFV, &GLOBAL::IONE);
		length = r * NN;
		dscal_(&length, &ntwo, EGFTU, &GLOBAL::IONE);
	};

	void EucBlindDeconvolution::EucHessianEta(Variable *x, Vector *etax, Vector *xix) const
	{
		const double *xptr = x->ObtainReadData();
		const double *Uptr = xptr;
		const double *Vptr = xptr + K * r;
		const double *dUptr = etax->ObtainReadData();
		const double *dVptr = dUptr + K * r;

		double *BU = new double[L * r * 4];
		double *CV = BU + L * r * 2;
		double *mdiagBetaxA = new double[L * 2];
		double *tmp = new double[4];

		for (integer i = 0; i < 2 * L; i++)
		{
			mdiagBetaxA[i] = 0;
		}

		integer LL2 = 2 * L, KK = K, NN = N, rr = r, LL = L;
		BtimesU(const_cast<double *> (dUptr), BU);
		fftwrapper(L, r, (fftw_complex *)BU, (fftw_complex *)BU, FFTW_FORWARD);
		CtimesV(const_cast<double *> (Vptr), CV);
		fftwrapper(L, r, (fftw_complex *)CV, (fftw_complex *)CV, FFTW_FORWARD);
		for (integer i = 0; i < L; i++)
		{
			dgemm_(GLOBAL::N, GLOBAL::T, &GLOBAL::ITWO, &GLOBAL::ITWO, &rr, &GLOBAL::DONE, BU + 2 * i, &LL2, CV + 2 * i, &LL2, &GLOBAL::DZERO, tmp, &GLOBAL::ITWO);
			mdiagBetaxA[2 * i] -= (tmp[0] - tmp[3]);
			mdiagBetaxA[2 * i + 1] -= (tmp[1] + tmp[2]);
		}

		BtimesU(const_cast<double *> (Uptr), BU);
		fftwrapper(L, r, (fftw_complex *)BU, (fftw_complex *)BU, FFTW_FORWARD);
		CtimesV(const_cast<double *> (dVptr), CV);
		fftwrapper(L, r, (fftw_complex *)CV, (fftw_complex *)CV, FFTW_FORWARD);
		for (integer i = 0; i < L; i++)
		{
			dgemm_(GLOBAL::N, GLOBAL::T, &GLOBAL::ITWO, &GLOBAL::ITWO, &rr, &GLOBAL::DONE, BU + 2 * i, &LL2, CV + 2 * i, &LL2, &GLOBAL::DZERO, tmp, &GLOBAL::ITWO);
			mdiagBetaxA[2 * i] -= (tmp[0] - tmp[3]);
			mdiagBetaxA[2 * i + 1] -= (tmp[1] + tmp[2]);
		}
		delete[] tmp;
		delete[] BU;

		double *xixptr = xix->ObtainWriteEntireData();
		double *EHFV = xixptr, *EHFTU = xixptr + K * r;
		double *tmpp = new double[2 * r * L];

		/*Computing EHFV*/
		CtimesV(const_cast<double *> (Vptr), tmpp);
		fftwrapper(L, r, (fftw_complex *)tmpp, (fftw_complex *)tmpp, FFTW_FORWARD);
		integer length = L * r;
		dscal_(&length, &GLOBAL::DNONE, tmpp + 1, &GLOBAL::ITWO);
#ifndef MATLAB_MEX_FILE
		for (integer i = 0; i < L; i++)
			zscal_(&rr, (doublecomplex *)(mdiagBetaxA + 2 * i), (doublecomplex *)(tmpp + 2 * i), &LL);
#else
		for (integer i = 0; i < L; i++)
			zscal_(&rr, mdiagBetaxA + 2 * i, tmpp + 2 * i, &LL);
#endif
		fftwrapper(L, r, (fftw_complex *)tmpp, (fftw_complex *)tmpp, FFTW_BACKWARD);
		BTtimesM(tmpp, EHFV);

		/*Computing EHFTU*/
		BtimesU(const_cast<double *> (Uptr), tmpp);
		fftwrapper(L, r, (fftw_complex *)tmpp, (fftw_complex *)tmpp, FFTW_FORWARD);
		length = L * r;
		dscal_(&length, &GLOBAL::DNONE, tmpp + 1, &GLOBAL::ITWO);
#ifndef MATLAB_MEX_FILE
		for (integer i = 0; i < L; i++)
			zscal_(&rr, (doublecomplex *)(mdiagBetaxA + 2 * i), (doublecomplex *)(tmpp + 2 * i), &LL);
#else
		for (integer i = 0; i < L; i++)
			zscal_(&rr, mdiagBetaxA + 2 * i, tmpp + 2 * i, &LL);
#endif
		fftwrapper(L, r, (fftw_complex *)tmpp, (fftw_complex *)tmpp, FFTW_BACKWARD);
		CTtimesM(tmpp, EHFTU);

		delete[] tmpp;
		delete[] mdiagBetaxA;

		length = r * KK;
		double ntwo = -2;
		dscal_(&length, &ntwo, EHFV, &GLOBAL::IONE);
		length = r * NN;
		dscal_(&length, &ntwo, EHFTU, &GLOBAL::IONE);

		const SharedSpace *ymdiagBXA = x->ObtainReadTempData("ymdiagBXA");
		const double *tmpptr = ymdiagBXA->ObtainReadData();
		double *EGFdV = new double[K * r], *EGFTdU = new double[N * r];
		tmpp = new double[2 * r * L];

		/*Computing EGFdV*/
		CtimesV(const_cast<double *> (dVptr), tmpp);
		fftwrapper(L, r, (fftw_complex *)tmpp, (fftw_complex *)tmpp, FFTW_FORWARD);
		length = L * r;
		dscal_(&length, &GLOBAL::DNONE, tmpp + 1, &GLOBAL::ITWO);
#ifndef MATLAB_MEX_FILE
		for (integer i = 0; i < L; i++)
			zscal_(&rr, (doublecomplex *)(tmpptr + 2 * i), (doublecomplex *)(tmpp + 2 * i), &LL);
#else
		for (integer i = 0; i < L; i++)
			zscal_(&rr, const_cast<double *> (tmpptr + 2 * i), tmpp + 2 * i, &LL);
#endif
		fftwrapper(L, r, (fftw_complex *)tmpp, (fftw_complex *)tmpp, FFTW_BACKWARD);
		BTtimesM(tmpp, EGFdV);

		/*Computing EGFTdU*/
		BtimesU(const_cast<double *> (dUptr), tmpp);
		fftwrapper(L, r, (fftw_complex *)tmpp, (fftw_complex *)tmpp, FFTW_FORWARD);
		length = L * r;
		dscal_(&length, &GLOBAL::DNONE, tmpp + 1, &GLOBAL::ITWO);
#ifndef MATLAB_MEX_FILE
		for (integer i = 0; i < L; i++)
			zscal_(&rr, (doublecomplex *)(tmpptr + 2 * i), (doublecomplex *)(tmpp + 2 * i), &LL);
#else
		for (integer i = 0; i < L; i++)
			zscal_(&rr, const_cast<double *> (tmpptr + 2 * i), tmpp + 2 * i, &LL);
#endif
		fftwrapper(L, r, (fftw_complex *)tmpp, (fftw_complex *)tmpp, FFTW_BACKWARD);
		CTtimesM(tmpp, EGFTdU);

		delete[] tmpp;
		ntwo = -2;
		length = r * KK;
		daxpy_(&length, &ntwo, EGFdV, &GLOBAL::IONE, EHFV, &GLOBAL::IONE);
		length = r * NN;
		daxpy_(&length, &ntwo, EGFTdU, &GLOBAL::IONE, EHFTU, &GLOBAL::IONE);

		delete[] EGFdV;
		delete[] EGFTdU;
	};

}; /*end of ROPTLIB namespace*/
#endif
