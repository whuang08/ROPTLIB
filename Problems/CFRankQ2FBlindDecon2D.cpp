
#include "Problems/CFRankQ2FBlindDecon2D.h"

#ifdef ROPTLIB_WITH_FFTW

/*Define the namespace*/
namespace ROPTLIB {

	CFRankQ2FBlindDecon2D::CFRankQ2FBlindDecon2D(Vector iny, SparseMatrix &inB, SparseMatrix &inC, integer inn1, integer inn2, integer inr, realdp inrho, realdp ind, realdp inmu)
	{
		y = iny;

        Bptr = &inB;
        Cptr = &inC;
        n1 = inn1;
        n2 = inn2;
        L = n1 * n2;
        r = inr;
        rho = inrho;
        d = ind;
        mu = inmu;
        
        NumGradHess = false;
        
		if (static_cast<int> (pow(2.0, static_cast<int>(log(static_cast<float> (L)) / log(static_cast<float> (2)) + 0.5))) - L != 0)
			printf("Warning: L must be a power of 2!\n");
	};

	CFRankQ2FBlindDecon2D::~CFRankQ2FBlindDecon2D(void)
	{
	};

	realdp CFRankQ2FBlindDecon2D::f(const Variable &x) const
	{
        Vector U = x.GetElement(0), V = x.GetElement(1);
        Vector FBU = ((*Bptr) * U).Reshape(n1, n2).GetFFT2D(FFTW_FORWARD).Reshape(L), FCV = ((*Cptr) * V).GetInvHaarFWT().Reshape(n1, n2).GetFFT2D(FFTW_BACKWARD).Reshape(L);
        const realdp *FBUptr = FBU.ObtainReadData(), *FCVptr = FCV.ObtainReadData();
        const realdp *yptr = y.ObtainReadData();
        Vector ymdiagBXC(L, "complex");
        realdp *ymdiagBXCptr = ymdiagBXC.ObtainWriteEntireData();
        
        for (integer i = 0; i < L; i++)
        {
            ymdiagBXCptr[2 * i] = yptr[2 * i] - (FBUptr[2 * i] * FCVptr[2 * i] + FBUptr[2 * i + 1] * FCVptr[2 * i + 1]);
            ymdiagBXCptr[2 * i + 1] = yptr[2 * i + 1] - (FBUptr[2 * i + 1] * FCVptr[2 * i] - FBUptr[2 * i] * FCVptr[2 * i + 1]);
        }

        realdp penalty = 0;
        if (rho != 0)
        {
            realdp tmp = static_cast<realdp> (L) / 8 / d / d / mu / mu;
            Vector rownorm2FBU(L + 1);
            realdp *rownorm2FBUptr = rownorm2FBU.ObtainWriteEntireData();
            integer length = r, LL2 = L * 2;
            for (integer i = 0; i < L; i++)
            {
                rownorm2FBUptr[i] = dot_(&length, const_cast<realdp *> (FBUptr) + i * 2, &LL2, const_cast<realdp *> (FBUptr) + i * 2, &LL2);
                rownorm2FBUptr[i] += dot_(&length, const_cast<realdp *> (FBUptr) + i * 2 + 1, &LL2, const_cast<realdp *> (FBUptr) + i * 2 + 1, &LL2);
            }
            length = L * 2 * r;
            const realdp *Vptr = V.ObtainReadData();
            rownorm2FBUptr[L] = dot_(&length, const_cast<realdp *> (Vptr), &GLOBAL::IONE, const_cast<realdp *> (Vptr), &GLOBAL::IONE);
            
            realdp tmp2 = 0;
            for (integer i = 0; i < L; i++)
            {
                tmp2 = ((tmp * rownorm2FBUptr[i] * rownorm2FBUptr[L] - 1 < 0) ? 0 : (tmp * rownorm2FBUptr[i] * rownorm2FBUptr[L] - 1));
                penalty += rho * tmp2 * tmp2;
            }
            x.AddToFields("rownorm2FBU", rownorm2FBU);
        }
        if(UseGrad)
        {
            x.AddToFields("ymdiagBXC", ymdiagBXC);
            x.AddToFields("FBU", FBU);
            x.AddToFields("FCV", FCV);
        }
        
        return ymdiagBXC.DotProduct(ymdiagBXC) + penalty;
	};

	Vector &CFRankQ2FBlindDecon2D::EucGrad(const Variable &x, Vector *result) const
	{
        result->NewMemoryOnWrite(); /*For multimanifolds, it is important to new memory at the beginning!*/
        
        Vector U = x.GetElement(0), V = x.GetElement(1);
        Vector ymdiagBXC = x.Field("ymdiagBXC"), FBU = x.Field("FBU"), FCV = x.Field("FCV");
        realdpcomplex ntwo = {-2.0, 0};
        
        Vector EGFV = ntwo * (ymdiagBXC.GetDiagTimesM(FCV, GLOBAL::L).Reshape(n1, n2).GetFFT2D(FFTW_BACKWARD).Reshape(L).GetTranspose() * (*Bptr)).GetTranspose();
        
        Vector EGFTU = ntwo * (ymdiagBXC.GetConj().GetDiagTimesM(FBU, GLOBAL::L).Reshape(n1, n2).GetFFT2D(FFTW_FORWARD).Reshape(L).GetHaarFWT().GetTranspose() * (*Cptr)).GetTranspose();
        
        if(rho != 0)
        {
            Vector rownorm2FBU = x.Field("rownorm2FBU");
            Vector DFBU = FBU;
            const realdp *rownorm2FBUptr = rownorm2FBU.ObtainReadData();
            realdp *DFBUptr = DFBU.ObtainWritePartialData();
            realdpcomplex coef = {0, 0};
            realdp tmp = static_cast<realdp> (L) / 8 / d / d / mu / mu;
            for (integer i = 0; i < L; i++)
            {
                coef.r = 2.0 * tmp * rho * 2 * ((tmp * rownorm2FBUptr[i] * rownorm2FBUptr[L] - 1 < 0) ? 0 : tmp * rownorm2FBUptr[i] * rownorm2FBUptr[L] - 1) * rownorm2FBUptr[L];
                scal_(&r, &coef, (realdpcomplex *)(DFBUptr + 2 * r * i), &L);
            }
            EGFV = EGFV + (DFBU.Reshape(n1, n2).GetFFT2D(FFTW_BACKWARD).Reshape(L).GetTranspose() * (*Bptr)).GetTranspose();
            
            realdpcomplex coef2 = {0, 0};
            for (integer i = 0; i < L; i++)
            {
                coef2.r += 2.0 * tmp * rho * rownorm2FBUptr[i] * 2 * ((tmp * rownorm2FBUptr[i] * rownorm2FBUptr[L] - 1 < 0) ? 0 : tmp * rownorm2FBUptr[i] * rownorm2FBUptr[L] - 1);
            }
            EGFTU = EGFTU + coef2 * V;
        }
        result->GetElement(0) = EGFV;
        result->GetElement(1) = EGFTU;

        return *result;
	};

	Vector &CFRankQ2FBlindDecon2D::EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const
	{ /*The Hessian of the penalty term is not considered here since the penalty is only C^1, not C^2 */
        result->NewMemoryOnWrite();/*For multimanifolds, it is important to new memory at the beginning!*/
        
        Vector U = x.GetElement(0), V = x.GetElement(1);
        Vector dU = etax.GetElement(0), dV = etax.GetElement(1);
        
        Vector ymdiagBXC = x.Field("ymdiagBXC"), FBU = x.Field("FBU"), FCV = x.Field("FCV");
        Vector FBdU = ((*Bptr) * dU).Reshape(n1, n2).GetFFT2D(FFTW_FORWARD).Reshape(L), FCdV = ((*Cptr) * dV).GetInvHaarFWT().Reshape(n1, n2).GetFFT2D(FFTW_BACKWARD).Reshape(L);
        realdpcomplex ntwo = {-2.0, 0}, two = {2.0, 0};
        
        Vector dDBUVC = FBU.GetHadamardProduct(FCdV.GetConj()) + FBdU.GetHadamardProduct(FCV.GetConj());
        
        Vector dJ1V = two * (dDBUVC.GetDiagTimesM(FCV, GLOBAL::L).Reshape(n1, n2).GetFFT2D(FFTW_BACKWARD).Reshape(L).GetTranspose() * (*Bptr)).GetTranspose();
        Vector J1dV = ntwo * (ymdiagBXC.GetDiagTimesM(FCdV, GLOBAL::L).Reshape(n1, n2).GetFFT2D(FFTW_BACKWARD).Reshape(L).GetTranspose() * (*Bptr)).GetTranspose();
        result->GetElement(0) = dJ1V + J1dV;
        
        Vector dJ1TU = two * (dDBUVC.GetConj().GetDiagTimesM(FBU, GLOBAL::L).Reshape(n1, n2).GetFFT2D(FFTW_FORWARD).Reshape(L).GetHaarFWT().GetTranspose() * (*Cptr)).GetTranspose();
        Vector J1TdU = ntwo * (ymdiagBXC.GetConj().GetDiagTimesM(FBdU, GLOBAL::L).Reshape(n1, n2).GetFFT2D(FFTW_FORWARD).Reshape(L).GetHaarFWT().GetTranspose() * (*Cptr)).GetTranspose();
        result->GetElement(1) = dJ1TU + J1TdU;
        
        return *result;
	};

}; /*end of ROPTLIB namespace*/
#endif
