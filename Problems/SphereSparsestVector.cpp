
#include "Problems/SphereSparsestVector.h"

/*Define the namespace*/
namespace ROPTLIB{

	SphereSparsestVector::SphereSparsestVector(Vector inQ)
	{
		Q = inQ;
		m = inQ.Getrow();
		n = inQ.Getcol();
        
        NumGradHess = false;
	};

	SphereSparsestVector::~SphereSparsestVector(void)
	{
	};

	realdp SphereSparsestVector::f(const Variable &x) const
	{
        Vector Qx(m); Qx.AlphaABaddBetaThis(1, Q, GLOBAL::N, x, GLOBAL::N, 0);
        Vector signQx(m);
        realdp *signQxptr = signQx.ObtainWriteEntireData();
        const realdp *Qxptr = Qx.ObtainReadData();
        
		realdp result = 0;
		for (integer i = 0; i < m; i++)
		{
			result += (Qxptr[i] < 0 ? -Qxptr[i] : Qxptr[i]);
			signQxptr[i] = static_cast<realdp> ((Qxptr[i] < 0 ? -1 : 1));
		}

		x.AddToFields("signQx", signQx);

		return result;
	};

	Vector &SphereSparsestVector::EucGrad(const Variable &x, Vector *result) const
	{
        Vector signQx = x.Field("signQx");
        result->AlphaABaddBetaThis(1, Q, GLOBAL::T, signQx, GLOBAL::N, 0);
        return *result;
	};
	
	Vector &SphereSparsestVector::EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const
	{
        *result = etax;
        return *result;
	};

}; /*end of ROPTLIB namespace*/
