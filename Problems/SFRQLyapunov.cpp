
#include "Problems/SFRQLyapunov.h"

/*Define the namespace*/
namespace ROPTLIB {

	SFRQLyapunov::SFRQLyapunov(SparseMatrix &inA, SparseMatrix &inM, Vector inC, integer inp)
	{
        Aptr = &inA;
        Mptr = &inM;
        C = inC;
        Cp = C.Getcol();
        n = C.Getrow();
        p = inp;
        NumGradHess = false;
	};

	SFRQLyapunov::~SFRQLyapunov(void)
	{
	};

	realdp SFRQLyapunov::f(const Variable &x) const
	{
        Vector AY = (*Aptr) * x;
        Vector MY = (*Mptr) * x;
        Vector CTY(Cp, p); CTY.AlphaABaddBetaThis(1, C, GLOBAL::T, x, GLOBAL::N, 0); /* CTY = C.GetTranspose() * x; */
        Vector YTAY(p, p); YTAY.AlphaABaddBetaThis(1, x, GLOBAL::T, AY, GLOBAL::N, 0); /* YTAY = x.GetTranspose() * AY; */
        Vector YTMY(p, p); YTMY.AlphaABaddBetaThis(1, x, GLOBAL::T, MY, GLOBAL::N, 0); /* YTMY = x.GetTranspose() * MY; */
        realdp result = YTAY.DotProduct(YTMY) - CTY.DotProduct(CTY);
        x.AddToFields("AY", AY);
        x.AddToFields("MY", MY);
        x.AddToFields("CTY", CTY);
        x.AddToFields("YTAY", YTAY);
        x.AddToFields("YTMY", YTMY);
        return result;
    };

	Vector &SFRQLyapunov::EucGrad(const Variable &x, Vector *result) const
	{
        Vector AY = x.Field("AY"), MY = x.Field("MY"), CTY = x.Field("CTY"), YTAY = x.Field("YTAY"), YTMY = x.Field("YTMY");
        result->AlphaABaddBetaThis(2, AY, GLOBAL::N, YTMY, GLOBAL::N, 0);
        result->AlphaABaddBetaThis(-2, C, GLOBAL::N, CTY, GLOBAL::N, 1);
        result->AlphaABaddBetaThis(2, MY, GLOBAL::N, YTAY, GLOBAL::N, 1); /* 2 * AY * YTMY - 2 * C * CTY + 2 * MY * YTAY; */
        return *result;
	};

	Vector &SFRQLyapunov::EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const
	{
        Vector AY = x.Field("AY"), MY = x.Field("MY"), YTAY = x.Field("YTAY"), YTMY = x.Field("YTMY");
        Vector Meta = (*Mptr) * etax;
        Vector YTMeta(p, p); YTMeta.AlphaABaddBetaThis(1, x, GLOBAL::T, Meta, GLOBAL::N, 0); /* YTMeta = x.GetTranspose() * Meta; */
        Vector Aeta = (*Aptr) * etax;
        Vector YTAeta(p, p); YTAeta.AlphaABaddBetaThis(1, x, GLOBAL::T, Aeta, GLOBAL::N, 0); /* YTAeta = x.GetTranspose() * Aeta; */
        
        result->AlphaABaddBetaThis(2, Aeta, GLOBAL::N, YTMY, GLOBAL::N, 0);
        result->AlphaABaddBetaThis(2, Meta, GLOBAL::N, YTAY, GLOBAL::N, 1);
        result->AlphaABaddBetaThis(2, AY, GLOBAL::N, YTMeta, GLOBAL::T, 1);
        result->AlphaABaddBetaThis(2, MY, GLOBAL::N, YTAeta, GLOBAL::T, 1);
        result->AlphaABaddBetaThis(2, AY, GLOBAL::N, YTMeta, GLOBAL::N, 1);
        result->AlphaABaddBetaThis(2, MY, GLOBAL::N, YTAeta, GLOBAL::N, 1);
        Vector tmp(Cp, p); tmp.AlphaABaddBetaThis(1, C, GLOBAL::T, etax, GLOBAL::N, 0);
        result->AlphaABaddBetaThis(-2, C, GLOBAL::N, tmp, GLOBAL::N, 1); /* 2 * (Aeta * YTMY + Meta * YTAY + AY * YTMeta.GetTranspose() + MY * YTAeta.GetTranspose() + AY * YTMeta + MY * YTAeta - C * (C.GetTranspose() * etax)); */
        return *result;
	};
 
	Vector &SFRQLyapunov::PreConditioner(const Variable &x, const Vector &eta, Vector *result) const
	{
        *result = eta;
        return *result; // no preconditioner!
		return LinearCG(x, eta, result);
	};

	Vector &SFRQLyapunov::ActionEH(const Variable &x, const Vector &intreta, Vector *result) const
	{
        Vector AY = x.Field("AY"), MY = x.Field("MY"), YTAY = x.Field("YTAY"), YTMY = x.Field("YTMY");
        Vector eta(n, p); Domain->ObtainExtr(x, intreta, &eta);
        Vector Meta = (*Mptr) * eta;
        Vector YTMeta(p, p); YTMeta.AlphaABaddBetaThis(1, x, GLOBAL::T, Meta, GLOBAL::N, 0); /* YTMeta = x.GetTranspose() * Meta; */
        Vector Aeta = (*Aptr) * eta;
        Vector YTAeta(p, p); YTAeta.AlphaABaddBetaThis(1, x, GLOBAL::T, Aeta, GLOBAL::N, 0); /* YTAeta = x.GetTranspose() * Aeta; */
        
        eta.AlphaABaddBetaThis(2, Aeta, GLOBAL::N, YTMY, GLOBAL::N, 0);
        eta.AlphaABaddBetaThis(2, Meta, GLOBAL::N, YTAY, GLOBAL::T, 1);
        eta.AlphaABaddBetaThis(2, AY, GLOBAL::N, YTMeta, GLOBAL::T, 1);
        eta.AlphaABaddBetaThis(2, MY, GLOBAL::N, YTAeta, GLOBAL::T, 1); /* 2 * (Aeta * YTMY + Meta * YTAY + AY * YTMeta.GetTranspose() + MY * YTAeta.GetTranspose()); */
        Domain->EucGradToGrad(x, eta, this, &eta);
        return Domain->ObtainIntr(x, eta, result);
	};

	Vector &SFRQLyapunov::LinearCG(const Variable x, const Vector &intreta, Vector *result) const
	{ /* solve A y - b = 0, where A is given by ActionnEH and b = intreta */
        integer length = intreta.Getlength();
        Vector y(intreta), Ay(intreta), r(intreta), p(intreta), Ap(intreta);
        ActionEH(x, y, &Ay);
        Domain->ScalarVectorAddVector(x, -1, intreta, Ay, &r); /*r = A y - b*/
        Domain->ScalarTimesVector(x, -1, r, &p); /*p = -r*/
        
		integer k = 0, maxiter = length;
		realdp alpha = 0, beta = 0;
		realdp err = Domain->Metric(x, r, r), newerr = 0;
		realdp err0 = err;
		while (err / err0 > 1e-3 && k < maxiter)
		{
			ActionEH(x, p, &Ap);
			alpha = err / Domain->Metric(x, p, Ap);
            Domain->ScalarVectorAddVector(x, alpha, p, y, &y); /*y = y + alpha * p*/
            Domain->ScalarVectorAddVector(x, alpha, Ap, r, &r); /*r = r + alpha * Ap*/
            newerr = Domain->Metric(x, r, r);
            beta = newerr / err;
            Domain->VectorLinearCombination(x, -1, r, beta, p, &p);
			err = newerr;
			k++;
		}
        if(err / err0 > 1e-3)
        {
            printf("warning: SFRQLyapunov::LinearCG does not find accurate enough solution!\n");
        }
        *result = y;
        return *result;
	};
}; /*end of ROPTLIB namespace*/
