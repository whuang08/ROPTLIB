
#include "Problems/GrassRQ.h"

/*Define the namespace*/
namespace ROPTLIB{

	GrassRQ::GrassRQ(Vector inB, integer inn, integer inp)
	{
		B = inB;
		n = inn;
		p = inp;
        
        NumGradHess = false;
	};

	GrassRQ::~GrassRQ(void)
	{
	};

	realdp GrassRQ::f(const Variable &x) const
	{
        Vector Bx(n, p);
        Bx.AlphaABaddBetaThis(1, B, GLOBAL::N, x, GLOBAL::N, 0);
        x.AddToFields("Bx", Bx);
        return (Bx.DotProduct(x) / 2);
	};

	Vector &GrassRQ::EucGrad(const Variable &x, Vector *result) const
	{
        *result = x.Field("Bx");
        return *result;
	};

	Vector &GrassRQ::EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const
	{
        result->AlphaABaddBetaThis(1, B, GLOBAL::N, etax, GLOBAL::N, 0);
        return *result;
	};
}; /*end of ROPTLIB namespace*/
