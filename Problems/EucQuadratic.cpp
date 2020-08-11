
#include "Problems/EucQuadratic.h"

/*Define the namespace*/
namespace ROPTLIB{

	EucQuadratic::EucQuadratic(Vector M)
	{
		A = M;
        
        NumGradHess = false;
	};

	EucQuadratic::~EucQuadratic(void)
	{
	};

	realdp EucQuadratic::f(const Variable &x) const
	{
        Vector Ax(x); Ax.AlphaABaddBetaThis(1, A, GLOBAL::N, x, GLOBAL::N, 0); /* Ax = A * x; */
        realdp result = 0.5 * (x.DotProduct(Ax));
        x.AddToFields("Ax", Ax);
        return result;
	};

	Vector &EucQuadratic::EucGrad(const Variable &x, Vector *result) const
	{
        *result = x.Field("Ax");
        return *result;
	};

	Vector &EucQuadratic::EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const
	{
        result->AlphaABaddBetaThis(1, A, GLOBAL::N, etax, GLOBAL::N, 0);
        return *result;
	};
}; /*end of ROPTLIB namespace*/
