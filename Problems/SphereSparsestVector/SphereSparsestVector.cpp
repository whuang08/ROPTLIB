
#include "Problems/SphereSparsestVector/SphereSparsestVector.h"

/*Define the namespace*/
namespace ROPTLIB{

	SphereSparsestVector::SphereSparsestVector(double *inQ, integer inm, integer inn)
	{
		Q = inQ;
		m = inm;
		n = inn;
	};

	SphereSparsestVector::~SphereSparsestVector(void)
	{
	};

	double SphereSparsestVector::f(Variable *x) const
	{
		const double *xptr = x->ObtainReadData();
		SharedSpace *SharedQx = new SharedSpace(1, m);
		double *Qx = SharedQx->ObtainWriteEntireData();

		dgemm_(GLOBAL::N, GLOBAL::N, &m, &GLOBAL::IONE, &n, &GLOBAL::DONE, Q, &m, const_cast<double *> (xptr), &n, &GLOBAL::DZERO, Qx, &m);

		double result = 0;
		for (integer i = 0; i < m; i++)
		{
			result += (Qx[i] < 0 ? -Qx[i] : Qx[i]);
			Qx[i] = (Qx[i] < 0 ? -1 : 1);
		}

		x->AddToTempData("signQx", SharedQx);

		return result;
	};

	void SphereSparsestVector::EucGrad(Variable *x, Vector *egf) const
	{
		const SharedSpace *SharedQx = x->ObtainReadTempData("signQx");
		const double *Qx = SharedQx->ObtainReadData();
		double *egfptr = egf->ObtainWriteEntireData();
		dgemm_(GLOBAL::T, GLOBAL::N, &n, &GLOBAL::IONE, &m, &GLOBAL::DONE, Q, &m, const_cast<double *> (Qx), &m, &GLOBAL::DZERO, egfptr, &n);
	};
	
	void SphereSparsestVector::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const
	{
		etax->CopyTo(exix);
	};

}; /*end of ROPTLIB namespace*/
