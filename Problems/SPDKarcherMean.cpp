
#include "Problems/SPDKarcherMean.h"

/*Define the namespace*/
namespace ROPTLIB{

	SPDKarcherMean::SPDKarcherMean(Vector inLs, integer inn, integer innum)
	{
		Ls = inLs;
		n = inn;
		num = innum;
        NumGradHess = false;
	};

	SPDKarcherMean::~SPDKarcherMean(void)
	{
	};

	realdp SPDKarcherMean::f(const Variable &x) const
	{
        if (!x.FieldsExist("_L"))
        {
            /* ideally x is a symmetric matrix and this tmp data is not needed.
            But when computing numerical gradient, the x may not be symmetric. This symmetrization is for
            the correctness of the numerical gradient*/
            Vector tmp = (x + x.GetTranspose()) / 2;
            tmp.CholDecom();
            x.AddToFields("_L", tmp.Field("_L"));
        }
        
        Vector SharedElement(n, n);
        Vector SharedLXL(1, &SharedElement, num);
        SharedLXL.NewMemoryOnWrite();
        Vector tmp1(n, n), tmp2(n, n);
        for(integer i = 0; i < num; i++)
        {
            tmp1 = x.Field("_L").TriangleLinSol(Ls.GetElement(i));
            tmp2.AlphaABaddBetaThis(1, tmp1, GLOBAL::N, tmp1, GLOBAL::T, 0);
            SharedLXL.GetElement(i) = tmp2.LogSym();
            SharedLXL.GetElement(i).AddToFields("_EigVec", tmp2.Field("_EigVec"));
            SharedLXL.GetElement(i).AddToFields("_EigVal", tmp2.Field("_EigVal"));
        }
        
        realdp result = SharedLXL.DotProduct(SharedLXL) / (static_cast<realdp> (2) * num);
        x.AddToFields("logLXL", SharedLXL);
        return result;
	};

	Vector &SPDKarcherMean::EucGrad(const Variable &x, Vector *result) const
	{
        Vector SharedLXL = x.Field("logLXL");
        result->SetToZeros();
        Vector tmp(n, n);
        for(integer i = 0; i < num; i++)
        {
            tmp.AlphaABaddBetaThis(1, SharedLXL.GetElement(i), GLOBAL::N, Ls.GetElement(i), GLOBAL::T, 0);
            result->AlphaXaddThis(1, tmp.TriangleLinSol(Ls.GetElement(i), GLOBAL::T));
        }
        
        *result = result->GetTranspose().TriangleLinSol(x.Field("_L"), GLOBAL::N).TriangleLinSol(x.Field("_L"), GLOBAL::T).GetTranspose();
        result->ScalarTimesThis(1.0 / num);
        return *result;
	};

	Vector &SPDKarcherMean::EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const
	{ /* DISCRETE REGRESSION METHODS ON THE CONE OF POSITIVE-DEFINITE MATRICES */
        result->SetToZeros();
        Vector EGrad = x.Field("EGrad");
        Vector SharedLXL = x.Field("logLXL");
        Vector U(n, n), Lambda(n), Dir(n, n), H(n, n), Z(n, n), tmp(n, n), tmp2(n, n);
        realdp *Zptr = Z.ObtainWriteEntireData();
        for(integer i = 0; i < num; i++)
        {
            tmp = etax.GetTranspose().TriangleLinSol(Ls.GetElement(i)).GetTranspose();
            Dir = tmp.TriangleLinSol(Ls.GetElement(i)); /* Dir = Ls.GetElement(i) % (etax / Ls.GetElement(i).GetTranspose()); */
            
            U = SharedLXL.GetElement(i).Field("_EigVec");
            Lambda = SharedLXL.GetElement(i).Field("_EigVal");
            
            tmp.AlphaABaddBetaThis(1, U, GLOBAL::T, Dir, GLOBAL::N, 0);
            H.AlphaABaddBetaThis(1, tmp, GLOBAL::N, U, GLOBAL::N, 0); /* H = U.GetTranspose() * Dir * U; */
            
            realdp *Lambdaptr = Lambda.ObtainWritePartialData();
            for(integer j = 0; j < n; j++)
                for(integer k = 0; k < n; k++)
                    Zptr[j + k * n] = (Lambdaptr[j] == Lambdaptr[k]) ? 1.0 / Lambdaptr[j] : ((std::log(Lambdaptr[j]) - std::log(Lambdaptr[k])) / (Lambdaptr[j] - Lambdaptr[k]));
            
            tmp.AlphaABaddBetaThis(1, H.GetHadamardProduct(Z), GLOBAL::N, U, GLOBAL::T, 0);
            tmp2.AlphaABaddBetaThis(1, U, GLOBAL::N, tmp, GLOBAL::N, 0); /* U * H.GetHadamardProduct(Z) * U.GetTranspose() */
            tmp = tmp2.TriangleLinSol(Ls.GetElement(i), GLOBAL::T);
            result->AlphaABaddBetaThis(1, tmp, GLOBAL::N, Ls.GetElement(i), GLOBAL::T, 1); /* *result = *result + Ls.GetElement(i).GetTranspose() % (U * H.GetHadamardProduct(Z) * U.GetTranspose()) * Ls.GetElement(i).GetTranspose(); */
        }
        *result = (*result / num - EGrad * etax) / x;
        return *result;
	};
}; /*end of ROPTLIB namespace*/
