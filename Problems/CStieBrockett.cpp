
#include "Problems/CStieBrockett.h"

/*Define the namespace*/
namespace ROPTLIB{

	CStieBrockett::CStieBrockett(Vector inB, Vector inD)
	{
        B = inB;
        D = inD;
        n = B.Getrow();
        p = D.Getrow();
        
        NumGradHess = false;
	};

	CStieBrockett::~CStieBrockett(void)
	{
	};

	realdp CStieBrockett::f(const Variable &x) const
	{
        Vector BxD(n, p, "complex"); BxD.AlphaABaddBetaThis(GLOBAL::ZONE, B, GLOBAL::N, x, GLOBAL::N, GLOBAL::ZZERO);
        D.DiagTimesM(BxD, GLOBAL::R); /* BxD = B * x; BxD = D.GetDiagTimesM(BxD, "R");*/
        realdp result = x.DotProduct(BxD);
        x.AddToFields("BxD", BxD);
        return result;
	};

	Vector &CStieBrockett::EucGrad(const Variable &x, Vector *result) const
	{
        *result = x.Field("BxD");
        result->ScalarTimesThis(2);
        return *result;
	};

	Vector &CStieBrockett::EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const
	{
        result->AlphaABaddBetaThis(GLOBAL::ZTWO, B, GLOBAL::N, etax, GLOBAL::N, GLOBAL::ZZERO);
        D.DiagTimesM(*result, GLOBAL::R);
        return *result; /* 2 * D.GetDiagTimesM(B * etax, "R"); */
	};
}; /*end of ROPTLIB namespace*/
