
#include "Problems/FRankEWeightApprox.h"

/*Define the namespace*/
namespace ROPTLIB{

	FRankEWeightApprox::FRankEWeightApprox(Vector inA, Vector inW, integer inm, integer inn, integer inr)
	{
		A = inA;
		W = inW;
		m = inm;
		n = inn;
		r = inr;
        
        NumGradHess = false;
	};

	FRankEWeightApprox::~FRankEWeightApprox(void)
	{
	};

	realdp FRankEWeightApprox::f(const Variable &x) const
	{
        Vector xmA(x); xmA.AlphaXaddThis(-1, A); /* xmA = x - A; */
        Vector WxmA(m * n, 1); WxmA.AlphaABaddBetaThis(1, W, GLOBAL::N, xmA.Reshape(m * n, 1), GLOBAL::N, 0);
        WxmA.Reshape(m, n); /* WxmA = (W * xmA.Reshape(m * n)).Reshape(m, n); */
        x.AddToFields("xmA", xmA);
        x.AddToFields("WxmA", WxmA);
        return xmA.DotProduct(WxmA);
	};
	
	Vector &FRankEWeightApprox::EucGrad(const Variable &x, Vector *result) const
	{
        *result = x.Field("WxmA"); result->ScalarTimesThis(2);
        return *result;
	};
	
	Vector &FRankEWeightApprox::EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const
	{
        Vector etaxreshape = etax; etaxreshape.Reshape(m * n, 1);
        result->Reshape(m * n, 1); result->AlphaABaddBetaThis(2.0, W, GLOBAL::N, etaxreshape, GLOBAL::N, 0);
        result->Reshape(m, n);
        return *result; /* 2.0 * (W * etax.Reshape(m * n)).Reshape(m, n); */
	};
}; /*end of ROPTLIB namespace*/
