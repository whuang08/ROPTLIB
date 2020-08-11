
#include "Problems/FRankESparseApprox.h"

/*Define the namespace*/
namespace ROPTLIB{

	FRankESparseApprox::FRankESparseApprox(Vector inA, realdp inlambda, integer inm, integer inn, integer inr, integer inlengthW)
	{
		A = inA;
        lambda = inlambda;
		m = inm;
		n = inn;
		r = inr;
        lengthW = inlengthW;
        L = 2;
        
        NumGradHess = false;
	};

	FRankESparseApprox::~FRankESparseApprox(void)
	{
	};

	realdp FRankESparseApprox::f(const Variable &x) const
	{
        Vector xmA(x); xmA.AlphaXaddThis(-1, A); /* xmA = x - A; */
        realdp result = xmA.DotProduct(xmA);
        
        const realdp *xptr = x.ObtainReadData();
        for(integer i = 0; i < n * m; i++)
            result += lambda * std::fabs(xptr[i]);
        
        x.AddToFields("xmA", xmA);
        return result;
	};
	
	Vector &FRankESparseApprox::EucGrad(const Variable &x, Vector *result) const
	{
        *result = x.Field("xmA"); result->ScalarTimesThis(2);
        return *result;
	};
	
	Vector &FRankESparseApprox::EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const
	{
        *result = etax; result->ScalarTimesThis(2);
        return *result;
	};

    Vector &FRankESparseApprox::PreConditioner(const Variable &x, const Vector &eta, Vector *result) const
    {
        Vector exresult(1);
        realdp *exresultptr = exresult.ObtainWriteEntireData();
        exresultptr[0] = L;
        *result = exresult;
        return *result;
    };

    Vector &FRankESparseApprox::ProxW(const Vector &x, const Vector &Weight, Vector *result) const
    { /* output = min(0, X + mu ./ W) + max(0, X - mu ./ W) */
        realdp *resultptr = result->ObtainWriteEntireData();
        const realdp *xptr = x.ObtainReadData();
        const realdp *Wptr = Weight.ObtainReadData();
        
        
        integer nblock = Weight.Getlength();
        integer blocksize = x.Getlength() / nblock;
        integer idx = 0;
        
        for(integer i = 0; i < nblock; i++)
        {
            for(integer j = 0; j < blocksize; j++)
            {
                idx = j + i * blocksize;
                resultptr[idx] = ((xptr[idx] + lambda / Wptr[i] < 0.0) ? xptr[idx] + lambda / Wptr[i] : 0.0) + ((xptr[idx] - lambda / Wptr[i] > 0.0) ? xptr[idx] - lambda / Wptr[i] : 0.0);
            }
        }
        return *result;
    };

    Vector &FRankESparseApprox::CalJW(const Vector &x, const Vector &eta, const Vector &Weight, Vector *result) const
    { /* output = (abs(x) .* W > lambda) .* eta; */
        const realdp *xptr = x.ObtainReadData();
        const realdp *etaptr = eta.ObtainReadData();
        realdp *resultptr = result->ObtainWriteEntireData();
        const realdp *Wptr = Weight.ObtainReadData();
        
        integer nblock = Weight.Getlength();
        integer blocksize = x.Getlength() / nblock;
        integer idx = 0;
        
        for(integer i = 0; i < nblock; i++)
        {
            for(integer j = 0; j < blocksize; j++)
            {
                idx = j + i * blocksize;
                resultptr[idx] = (std::abs(xptr[idx]) * Wptr[i] > lambda) ? etaptr[idx] : 0;
            }
        }
        return *result;
    };
}; /*end of ROPTLIB namespace*/
