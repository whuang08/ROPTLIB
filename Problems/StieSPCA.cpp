
#include "Problems/StieSPCA.h"

/*Define the namespace*/
namespace ROPTLIB{

	StieSPCA::StieSPCA(Vector inB, realdp inlambda, integer inn, integer inm, integer inp, integer inlengthW)
	{
		B = inB;
		n = inn;
        m = inm;
		p = inp;
        lengthW = inlengthW;
        lambda = inlambda;
        colnormB = Vector (n);
        
        realdp *colnormBptr = colnormB.ObtainWriteEntireData();
        const realdp *Bptr = B.ObtainReadData();
        
        for(integer i = 0; i < n; i++)
        {
            colnormBptr[i] = dot_(&m, const_cast<realdp *> (Bptr + i * m), &GLOBAL::IONE, const_cast<realdp *> (Bptr + i * m), &GLOBAL::IONE);
        }
        
        B.SVDDecom();
        Vector S = B.Field("_S");
        L = S.ObtainReadData()[0] * S.ObtainReadData()[0] * 2;
        B.RemoveAllFromFields();
        
        NumGradHess = false;
	};

	StieSPCA::~StieSPCA(void)
	{
	};

	realdp StieSPCA::f(const Variable &x) const
	{
        Vector Bx(m, p); Bx.AlphaABaddBetaThis(1, B, GLOBAL::N, x, GLOBAL::N, 0);/* Bx = B * x; */
        realdp result = - Bx.DotProduct(Bx);
        
		const realdp *xptr = x.ObtainReadData();
        for(integer i = 0; i < n * p; i++)
            result += lambda * std::fabs(xptr[i]);
        x.AddToFields("Bx", Bx);
        return result;
	};

	Vector &StieSPCA::EucGrad(const Variable &x, Vector *result) const
	{
        Vector Bx = x.Field("Bx");
        result->AlphaABaddBetaThis(-2, B, GLOBAL::T, Bx, GLOBAL::N, 0); /* -2.0 * (B.GetTranspose() * Bx) */
        
        return *result;
	};

	Vector &StieSPCA::EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const
	{
        Vector Betax(m, p); Betax.AlphaABaddBetaThis(1, B, GLOBAL::N, etax, GLOBAL::N, 0);
        result->AlphaABaddBetaThis(-2, B, GLOBAL::T, Betax, GLOBAL::N, 0); /* -2.0 * (B.GetTranspose() * (B * etax)) */
        
        return *result;
	};

    Vector &StieSPCA::PreConditioner(const Variable &x, const Vector &eta, Vector *result) const
    {
        if(lengthW == x.Getlength())
        {
            Vector Bx = x.Field("Bx");
            const realdp *Bxptr = Bx.ObtainReadData();
            const realdp *colnormBptr = colnormB.ObtainReadData();
            realdp *w1 = new realdp[p];
            for (integer i = 0; i < p; i++)
            {
                w1[i] = dot_(&m, const_cast<realdp *>(Bxptr + i * m), &GLOBAL::IONE, const_cast<realdp *>(Bxptr + i * m), &GLOBAL::IONE);
            }
            *result = Vector (n, p);
            realdp *resultptr = result->ObtainWriteEntireData();
            for(integer i = 0; i < n; i++)
                for(integer j = 0; j < p; j++)
                    resultptr[i + j * n] = ((w1[j] - colnormBptr[i]) * 2.0 > 0.1) ? (w1[j] - colnormBptr[i]) * 2.0 : 0.1;
            delete [] w1;
            return *result;
        } else
        if(lengthW == p)
        {
            Vector Bx = x.Field("Bx");
            const realdp *Bxptr = Bx.ObtainReadData();
            *result = Vector (p);
            realdp *resultptr = result->ObtainWriteEntireData();
            for (integer i = 0; i < p; i++)
            {
                resultptr[i] = dot_(&m, const_cast<realdp *>(Bxptr + i * m), &GLOBAL::IONE, const_cast<realdp *>(Bxptr + i * m), &GLOBAL::IONE);
                resultptr[i] = ((resultptr[i] - 1.0) * 2.0 > 0.1) ? (resultptr[i] - 1.0) * 2.0 : 0.1;
            }
            return *result;
        } else
        if(lengthW == 1)
        {
            *result = Vector (1);
            realdp *resultptr = result->ObtainWriteEntireData();
            resultptr[0] = L;
            return *result;
        } else
        {
            printf("Warning: StieSPCA::PreConditioner: the size of weighting matrix is not supported.\n");
            *result = Vector (1);
            result->SetToZeros();
            return *result;
        }
    };

    Vector &StieSPCA::ProxW(const Vector &x, const Vector &Weight, Vector *result) const
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

    Vector &StieSPCA::CalJW(const Vector &x, const Vector &eta, const Vector &Weight, Vector *result) const
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
