
#include "Problems/ProdStieSumBrockett.h"

/*Define the namespace*/
namespace ROPTLIB{

	ProdStieSumBrockett::ProdStieSumBrockett(Vector inB1, Vector inD1, Vector inB2, Vector inD2, Vector inB3, Vector inD3)
	{
		B1 = inB1;
		D1 = inD1;
		B2 = inB2;
		D2 = inD2;
		B3 = inB3;
		D3 = inD3;
        
        n = B1.Getrow();
        p = D1.Getlength();
        m = B3.Getrow();
        q = D3.Getlength();
        
        NumGradHess = false;
	};

	ProdStieSumBrockett::~ProdStieSumBrockett(void)
	{
	};

	realdp ProdStieSumBrockett::f(const Variable &x) const
	{
        Vector B1x1D1(n, p); B1x1D1.AlphaABaddBetaThis(1, B1, GLOBAL::N, x.GetElement(0), GLOBAL::N, 0); /* B1x1D1 = B1 * x.GetElement(0); */
        D1.DiagTimesM(B1x1D1, GLOBAL::R);
        realdp result = x.GetElement(0).DotProduct(B1x1D1);
        x.AddToFields("B1x1D1", B1x1D1);
        
        Vector B2x2D2(n, p); B2x2D2.AlphaABaddBetaThis(1, B2, GLOBAL::N, x.GetElement(1), GLOBAL::N, 0); /* B2x2D2 = B2 * x.GetElement(1); */
        D2.DiagTimesM(B2x2D2, GLOBAL::R);
        result += x.GetElement(1).DotProduct(B2x2D2);
        x.AddToFields("B2x2D2", B2x2D2);
        
        Vector B3x3D3(m, q); B3x3D3.AlphaABaddBetaThis(1, B3, GLOBAL::N, x.GetElement(2), GLOBAL::N, 0); /* B3x3D3 = B3 * x.GetElement(2); */
        D3.DiagTimesM(B3x3D3, GLOBAL::R);
        result += x.GetElement(2).DotProduct(B3x3D3);
        x.AddToFields("B3x3D3", B3x3D3);
        
        return result;
	};

	Vector &ProdStieSumBrockett::EucGrad(const Variable &x, Vector *result) const
	{
        result->NewMemoryOnWrite();
        result->GetElement(0) = x.Field("B1x1D1");
        result->GetElement(1) = x.Field("B2x2D2");
        result->GetElement(2) = x.Field("B3x3D3");
        Domain->ScalarTimesVector(x, 2, *result, result);
        return *result;
	};

	Vector &ProdStieSumBrockett::EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const
	{
        result->NewMemoryOnWrite();
        result->GetElement(0).AlphaABaddBetaThis(2, B1, GLOBAL::N, etax.GetElement(0), GLOBAL::N, 0);
        D1.DiagTimesM(result->GetElement(0), GLOBAL::R);
        result->GetElement(1).AlphaABaddBetaThis(2, B2, GLOBAL::N, etax.GetElement(1), GLOBAL::N, 0);
        D2.DiagTimesM(result->GetElement(1), GLOBAL::R);
        result->GetElement(2).AlphaABaddBetaThis(2, B3, GLOBAL::N, etax.GetElement(2), GLOBAL::N, 0);
        D3.DiagTimesM(result->GetElement(2), GLOBAL::R);
        
        return *result;
	};
}; /*end of ROPTLIB namespace*/
