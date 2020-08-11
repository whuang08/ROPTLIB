
#include "Problems/StieBrockett.h"

/*Define the namespace*/
namespace ROPTLIB{

	StieBrockett::StieBrockett(Vector inB, Vector inD)
	{
        B = inB;
        D = inD;
        n = B.Getrow();
        p = D.Getlength();
        
        NumGradHess = false;
	};

	StieBrockett::~StieBrockett(void)
	{
	};

	realdp StieBrockett::f(const Variable &x) const
	{
        Vector BxD(n, p); BxD.AlphaABaddBetaThis(1, B, GLOBAL::N, x, GLOBAL::N, 0);
        D.DiagTimesM(BxD, GLOBAL::R); /* BxD = B * x; BxD = D.GetDiagTimesM(BxD, "R");*/
        realdp result = x.DotProduct(BxD);
        x.AddToFields("BxD", BxD);
        return result;
	};

	Vector &StieBrockett::EucGrad(const Variable &x, Vector *result) const
	{
        *result = x.Field("BxD");
        result->ScalarTimesThis(2);
        return *result;
	};

	Vector &StieBrockett::EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const
	{
        result->AlphaABaddBetaThis(2, B, GLOBAL::N, etax, GLOBAL::N, 0);
        D.DiagTimesM(*result, GLOBAL::R);
        return *result; /* 2 * D.GetDiagTimesM(B * etax, "R"); */
	};
}; /*end of ROPTLIB namespace*/
