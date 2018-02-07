
#include "Problems/CSOPhaseRetrieval/CSOPhaseRetrieval.h"

#ifdef ROPTLIB_WITH_FFTW

/*Define the namespace*/
namespace ROPTLIB {

	CSOPhaseRetrieval::CSOPhaseRetrieval(double *inb, double *inmasks, double inkappa, integer inn1, integer inn2, integer inl, integer inr)
	{
		b = inb;
		masks = inmasks;
		kappa = inkappa;
		n1 = inn1;
		n2 = inn2;
		l = inl;
		r = inr;
		m = n1 * n2 * l;

		flags = FFTW_ESTIMATE;
		p = fftw_plan_dft_2d(n2, n1, nullptr, nullptr, FFTW_FORWARD, flags);
	};

	CSOPhaseRetrieval::~CSOPhaseRetrieval(void)
	{
		fftw_destroy_plan(p);
	};

	void CSOPhaseRetrieval::fftwrapper(integer inn1, integer inn2, fftw_complex *in, fftw_complex *out, int sign) const
	{
		p = fftw_plan_dft_2d(inn2, inn1, in, out, sign, flags);
		fftw_execute(p);
	};

	double CSOPhaseRetrieval::f(Variable *x) const
	{/*x is a 2 * n1 * n2 by r matrix*/
		SharedSpace *SharedcD = new SharedSpace(1, m);
		SharedSpace *SharedZY = new SharedSpace(2, 2 * m, r);
		double *cD = SharedcD->ObtainWriteEntireData();
		double *ZY = SharedZY->ObtainWriteEntireData();
		double sqn = sqrt(static_cast<double> (n1 * n2));
		
		for (integer i = 0; i < m; i++)
			cD[i] = 0;

		const double *xptr = x->ObtainReadData();
		double *sZqi = new double[n1 * n2];
		
		for (integer i = 0; i < l; i++)
		{
			for (integer j = 0; j < n1 * n2; j++)
				sZqi[j] = 0;

			for (integer j = 0; j < r; j++)
			{
				for (integer k = 0; k < n1 * n2; k++)
				{
					ZY[i * 2 * n1 * n2 + 2 * k + j * 2 * m] = (xptr[2 * k + j * 2 * n1 * n2] * masks[2 * k + i * 2 * n1 * n2] - xptr[2 * k + 1 + j * 2 * n1 * n2] * masks[2 * k + 1 + i * 2 * n1 * n2]) / sqn;
					ZY[i * 2 * n1 * n2 + 2 * k + 1 + j * 2 * m] = (xptr[2 * k + 1 + j * 2 * n1 * n2] * masks[2 * k + i * 2 * n1 * n2] + xptr[2 * k + j * 2 * n1 * n2] * masks[2 * k + 1 + i * 2 * n1 * n2]) / sqn;
				}
				fftwrapper(n1, n2, (fftw_complex *) (ZY + i * 2 * n1 * n2 + j * 2 * m), (fftw_complex *) (ZY + i * 2 * n1 * n2 + j * 2 * m), FFTW_FORWARD);
				//ForDebug::Print("ZY:", ZY, 2 * n1, n2);//--
				//std::cin >> cD[0];
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

		integer mm = m, nn2r = 2 * n1 * n2 * r;
		double result = ddot_(&mm, cD, &GLOBAL::IONE, cD, &GLOBAL::IONE) / ddot_(&mm, b, &GLOBAL::IONE, b, &GLOBAL::IONE);

		result = (kappa == 0) ? result : result + kappa * ddot_(&nn2r, const_cast<double *> (xptr), &GLOBAL::IONE, const_cast<double *> (xptr), &GLOBAL::IONE);

		x->AddToTempData("cD", SharedcD);
		x->AddToTempData("ZY", SharedZY);

		return result;
	};

	void CSOPhaseRetrieval::EucGrad(Variable *x, Vector *gf) const
	{
		const SharedSpace *SharedcD = x->ObtainReadTempData("cD");
		const SharedSpace *SharedZY = x->ObtainReadTempData("ZY");
		const double *cD = SharedcD->ObtainReadData();
		const double *ZY = SharedZY->ObtainReadData();
		double *gfptr = gf->ObtainWriteEntireData();
		for (integer i = 0; i < 2 * n1 * n2 * r; i++)
			gfptr[i] = 0;

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

		integer length = 2 * n1 * n2 * r;
		double sqn = sqrt(static_cast<double> (n1 * n2));
		for (integer i = 0; i < l; i++)
		{
			for (integer j = 0; j < r; j++)
			{
				fftwrapper(n1, n2, (fftw_complex *) (DZY + i * 2 * n1 * n2 + j * 2 * m), (fftw_complex *) (DZY + i * 2 * n1 * n2 + j * 2 * m), FFTW_BACKWARD);
				//ForDebug::Print("DZY:", DZY, 2 * n1, n2);//---
				for (integer k = 0; k < n1 * n2; k++)
				{
					double realv = (DZY[2 * k + i * 2 * n1 * n2 + j * 2 * m] * masks[2 * k + i * 2 * n1 * n2] + DZY[2 * k + 1 + i * 2 * n1 * n2 + j * 2 * m] * masks[2 * k + 1 + i * 2 * n1 * n2]) / sqn;
					double imagv = (-DZY[2 * k + i * 2 * n1 * n2 + j * 2 * m] * masks[2 * k + 1 + i * 2 * n1 * n2] + DZY[2 * k + 1 + i * 2 * n1 * n2 + j * 2 * m] * masks[2 * k + i * 2 * n1 * n2]) / sqn;
					tmp[2 * k + j * n1 * n2 * 2] = realv;
					tmp[2 * k + 1 + j * n1 * n2 * 2] = imagv;
				}
				//ForDebug::Print("tmp:", tmp, 2 * n1, n2);//---
				//std::cin >> DZY[0];//---
			}
			daxpy_(&length, &GLOBAL::DONE, tmp, &GLOBAL::IONE, gfptr, &GLOBAL::IONE);
		}
		delete[] DZY;

		integer mm = m;
		double scalar = 4.0 / ddot_(&mm, b, &GLOBAL::IONE, b, &GLOBAL::IONE);
		dscal_(&length, &scalar, gfptr, &GLOBAL::IONE);

		scalar = 2.0 * kappa;
		const double *xptr = x->ObtainReadData();
		if(kappa != 0)
			daxpy_(&length, &scalar, const_cast<double *>(xptr), &GLOBAL::IONE, gfptr, &GLOBAL::IONE);

		//TODO Y / (x'*x) use compute HHR in CpxNstQorth
		CpxNStQOrth::ComputeHHR(x);
		const SharedSpace *HHR = x->ObtainReadTempData("HHR");
		const double *ptrHHR = HHR->ObtainReadData();
		
		double *YH = new double[2 * n1 * n2 * r];

		/*YH = Y^H*/
		for (integer i = 0; i < n1 * n2; i++)
		{
			for (integer j = 0; j < r; j++)
			{
				YH[j * 2 + i * 2 * r] = gfptr[2 * i + j * 2 * n1 * n2];
				YH[j * 2 + 1 + i * 2 * r] = -gfptr[2 * i + 1 + j * 2 * n1 * n2];
			}
		}

		integer rr = r, nn = n1 * n2, info;
#ifndef MATLAB_MEX_FILE
		ztrtrs_(GLOBAL::U, GLOBAL::C, GLOBAL::N, &rr, &nn, (doublecomplex *)ptrHHR, &nn, (doublecomplex *)YH, &rr, &info);
		ztrtrs_(GLOBAL::U, GLOBAL::N, GLOBAL::N, &rr, &nn, (doublecomplex *)ptrHHR, &nn, (doublecomplex *)YH, &rr, &info);
#else
		ztrtrs_(GLOBAL::U, GLOBAL::C, GLOBAL::N, &rr, &nn, const_cast<double *> (ptrHHR), &nn, YH, &rr, &info);
		ztrtrs_(GLOBAL::U, GLOBAL::N, GLOBAL::N, &rr, &nn, const_cast<double *> (ptrHHR), &nn, YH, &rr, &info);
#endif
		/*YH = Y^H*/
		for (integer i = 0; i < n1 * n2; i++)
		{
			for (integer j = 0; j < r; j++)
			{
				gfptr[2 * i + j * 2 * n1 * n2] = YH[j * 2 + i * 2 * r];
				gfptr[2 * i + 1 + j * 2 * n1 * n2] = -YH[j * 2 + 1 + i * 2 * r];
			}
		}

		delete[] YH;
	};

}; /*end of ROPTLIB namespace*/
#endif
