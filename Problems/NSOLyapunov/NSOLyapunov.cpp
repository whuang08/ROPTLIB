
#include "Problems/NSOLyapunov/NSOLyapunov.h"

/*Define the namespace*/
namespace ROPTLIB {

	NSOLyapunov::NSOLyapunov(double *inA, double *inM, double *inC, integer inn, integer inp)
	{
		A = inA;
		M = inM;
		C = inC;
		n = inn;
		p = inp;
	};

	NSOLyapunov::~NSOLyapunov(void)
	{
	};

	double NSOLyapunov::f(Variable *x) const
	{
		const double *Y = x->ObtainReadData();
		SharedSpace *SharedAY = new SharedSpace(2, n, p);
		SharedSpace *SharedMY = new SharedSpace(2, n, p);
		SharedSpace *SharedCY = new SharedSpace(2, n, p);
		SharedSpace *SharedYTAY = new SharedSpace(2, p, p);
		SharedSpace *SharedYTMY = new SharedSpace(2, p, p);

		double *AY = SharedAY->ObtainWriteEntireData();
		double *MY = SharedMY->ObtainWriteEntireData();
		double *CY = SharedCY->ObtainWriteEntireData();
		double *YTAY = SharedYTAY->ObtainWriteEntireData();
		double *YTMY = SharedYTMY->ObtainWriteEntireData();

		integer N = n, P = p;
		/*AY = A * Y*/
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &N, &GLOBAL::DONE, A, &N, const_cast<double *>(Y), &N, &GLOBAL::DZERO, AY, &N);

		/*MY = M * Y*/
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &N, &GLOBAL::DONE, M, &N, const_cast<double *>(Y), &N, &GLOBAL::DZERO, MY, &N);

		/*CY = C * Y*/
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &N, &GLOBAL::DONE, C, &N, const_cast<double *>(Y), &N, &GLOBAL::DZERO, CY, &N);

		/*YTAY = Y^T * A * Y*/
		dgemm_(GLOBAL::T, GLOBAL::N, &P, &P, &N, &GLOBAL::DONE, const_cast<double *>(Y), &N, AY, &N, &GLOBAL::DZERO, YTAY, &P);

		/*YTMY = Y^T * M * Y*/
		dgemm_(GLOBAL::T, GLOBAL::N, &P, &P, &N, &GLOBAL::DONE, const_cast<double *>(Y), &N, MY, &N, &GLOBAL::DZERO, YTMY, &P);

		integer length = p * p, length2 = n * p;

		double result = ddot_(&length, YTAY, &GLOBAL::IONE, YTMY, &GLOBAL::IONE) - ddot_(&length2, const_cast<double *> (Y), &GLOBAL::IONE, CY, &GLOBAL::IONE);

		x->AddToTempData("AY", SharedAY);
		x->AddToTempData("MY", SharedMY);
		x->AddToTempData("CY", SharedCY);
		x->AddToTempData("YTAY", SharedYTAY);
		x->AddToTempData("YTMY", SharedYTMY);

		return result;
	};

	void NSOLyapunov::EucGrad(Variable *x, Vector *gf) const
	{
		const SharedSpace *SharedAY = x->ObtainReadTempData("AY");
		const SharedSpace *SharedMY = x->ObtainReadTempData("MY");
		const SharedSpace *SharedCY = x->ObtainReadTempData("CY");
		const SharedSpace *SharedYTAY = x->ObtainReadTempData("YTAY");
		const SharedSpace *SharedYTMY = x->ObtainReadTempData("YTMY");

		const double *AY = SharedAY->ObtainReadData();
		const double *MY = SharedMY->ObtainReadData();
		const double *CY = SharedCY->ObtainReadData();
		const double *YTAY = SharedYTAY->ObtainReadData();
		const double *YTMY = SharedYTMY->ObtainReadData();


		double *gfv = gf->ObtainWriteEntireData();

		integer length = n * p;
		dcopy_(&length, const_cast<double *>(CY), &GLOBAL::IONE, gfv, &GLOBAL::IONE);

		integer N = n, P = p;
		double NTWO = -2;
		/*gfv = 2 * AY * YTMY - 2 * C * Y*/
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &P, &GLOBAL::DTWO, const_cast<double *> (AY), &N, const_cast<double *> (YTMY), &P, &NTWO, gfv, &N);

		/*gfv = gfv + 2 * MY * YTAY*/
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &P, &GLOBAL::DTWO, const_cast<double *>(MY), &N, const_cast<double *>(YTAY), &P, &GLOBAL::DONE, gfv, &N);
	};

	void NSOLyapunov::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const
	{
		const SharedSpace *SharedAY = x->ObtainReadTempData("AY");
		const SharedSpace *SharedMY = x->ObtainReadTempData("MY");
		const SharedSpace *SharedYTAY = x->ObtainReadTempData("YTAY");
		const SharedSpace *SharedYTMY = x->ObtainReadTempData("YTMY");

		const double *Y = x->ObtainReadData();
		const double *AY = SharedAY->ObtainReadData();
		const double *MY = SharedMY->ObtainReadData();
		const double *YTAY = SharedYTAY->ObtainReadData();
		const double *YTMY = SharedYTMY->ObtainReadData();

		const double *etaTV = etax->ObtainReadData();
		double *exixTV = exix->ObtainWriteEntireData();

		double *Aeta = new double[2 * n * p + p * p * 2];
		double *Meta = Aeta + n * p;
		double *YTMeta = Meta + n * p;
		double *YTAeta = YTMeta + p * p;
		integer N = n, P = p;

		/*Meta = M * eta*/
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &N, &GLOBAL::DONE, M, &N, const_cast<double *> (etaTV), &N, &GLOBAL::DZERO, Meta, &N);
		/*YTMeta = Y^T * M * eta*/
		dgemm_(GLOBAL::T, GLOBAL::N, &P, &P, &N, &GLOBAL::DONE, const_cast<double *> (Y), &N, Meta, &N, &GLOBAL::DZERO, YTMeta, &P);
		/*Aeta = A * eta*/
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &N, &GLOBAL::DONE, A, &N, const_cast<double *> (etaTV), &N, &GLOBAL::DZERO, Aeta, &N);
		/*YTAeta = Y^T * M * eta*/
		dgemm_(GLOBAL::T, GLOBAL::N, &P, &P, &N, &GLOBAL::DONE, const_cast<double *> (Y), &N, Aeta, &N, &GLOBAL::DZERO, YTAeta, &P);

		/*exixTV = 2 * A * eta * Y^T * M * Y */
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &P, &GLOBAL::DTWO, Aeta, &N, const_cast<double *> (YTMY), &P, &GLOBAL::DZERO, exixTV, &N);

		/*exixTV = 2 * A * eta * Y^T * M * Y + 2 * M * eta * Y^T * A * Y*/
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &P, &GLOBAL::DTWO, Meta, &N, const_cast<double *> (YTAY), &P, &GLOBAL::DONE, exixTV, &N);

		/*exixTV = 2 * A * eta * Y^T * M * Y + 2 * M * eta * Y^T * A * Y + 2 * A * Y * eta^T * M * Y*/
		dgemm_(GLOBAL::N, GLOBAL::T, &N, &P, &P, &GLOBAL::DTWO, const_cast<double *> (AY), &N, YTMeta, &P, &GLOBAL::DONE, exixTV, &N);

		/*exixTV = 2 * A * eta * Y^T * M * Y + 2 * M * eta * Y^T * A * Y + 2 * A * Y * eta^T * M * Y + 2 * M * Y * eta^T * A * Y*/
		dgemm_(GLOBAL::N, GLOBAL::T, &N, &P, &P, &GLOBAL::DTWO, const_cast<double *> (MY), &N, YTAeta, &P, &GLOBAL::DONE, exixTV, &N);

		/*exixTV = 2 * A * eta * Y^T * M * Y + 2 * M * eta * Y^T * A * Y + 2 * A * Y * eta^T * M * Y + 2 * M * Y * eta^T * A * Y + 2 * A * Y * Y^T * M * \eta*/
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &P, &GLOBAL::DTWO, const_cast<double *> (AY), &N, YTMeta, &P, &GLOBAL::DONE, exixTV, &N);

		/*exixTV = 2 * A * eta * Y^T * M * Y + 2 * M * eta * Y^T * A * Y + 2 * A * Y * eta^T * M * Y + 2 * M * Y * eta^T * A * Y + 2 * A * Y * Y^T * M * \eta + 2 * M * Y * Y^T * A * eta*/
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &P, &GLOBAL::DTWO, const_cast<double *> (MY), &N, YTAeta, &P, &GLOBAL::DONE, exixTV, &N);

		/*exixTV = 2 * A * eta * Y^T * M * Y + 2 * M * eta * Y^T * A * Y + 2 * A * Y * eta^T * M * Y + 2 * M * Y * eta^T * A * Y + 2 * A * Y * Y^T * M * \eta + 2 * M * Y * Y^T * A * eta - 2 * C * \eta*/
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &N, &GLOBAL::DNTWO, C, &N, const_cast<double *> (etaTV), &N, &GLOBAL::DONE, exixTV, &N);

		delete[] Aeta;
	};
}; /*end of ROPTLIB namespace*/
