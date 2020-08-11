
#include "Problems/CSFRQPhaseRetrieval.h"

#ifdef ROPTLIB_WITH_FFTW

/*Define the namespace*/
namespace ROPTLIB {
	CSFRQPhaseRetrieval::CSFRQPhaseRetrieval(Vector inb, Vector inmasks, realdp inkappa, integer inn1, integer inn2, integer inl, integer inr)
	{
		b = inb;
		masks = inmasks;
		kappa = inkappa;
		n1 = inn1;
		n2 = inn2;
        n = n1 * n2;
		l = inl;
		r = inr;
		m = n1 * n2 * l;
        NumGradHess = false;
	};

	CSFRQPhaseRetrieval::~CSFRQPhaseRetrieval(void)
	{
	};

	realdp CSFRQPhaseRetrieval::f(const Variable &x) const
	{/*x is a n1 * n2 by r complex matrix*/
        Vector cD(m); cD.SetToZeros();
        Vector ZY(m, r, "complex");
        realdp sqn = sqrt(static_cast<realdp> (n1 * n2));
        Vector tmpv(n, "complex"), tmpv2(n, "complex");
        realdpcomplex *tmpvptr = (realdpcomplex *) tmpv.ObtainWriteEntireData();
        realdpcomplex *tmpv2ptr = (realdpcomplex *) tmpv2.ObtainWriteEntireData();
        realdpcomplex *ZYptr = (realdpcomplex *) ZY.ObtainWriteEntireData();
        const realdp *bptr = b.ObtainReadData();
        realdp *cDptr = cD.ObtainWriteEntireData();
        const realdpcomplex *masksptr = (realdpcomplex *) masks.ObtainReadData();
        const realdpcomplex *xptr = (realdpcomplex *) x.ObtainReadData();

        for(integer i = 0; i < l; i++)
        {
            for(integer j = 0; j < r; j++)
            {
                for(integer k = 0; k < n; k++)
                {
                    tmpvptr[k].r = (masksptr[k + i * n].r * xptr[k + j * n].r - masksptr[k + i * n].i * xptr[k + j * n].i) / sqn;
                    tmpvptr[k].i = (masksptr[k + i * n].i * xptr[k + j * n].r + masksptr[k + i * n].r * xptr[k + j * n].i) / sqn;
                }
                
                tmpv.Reshape(n1, n2).FFT2D(FFTW_FORWARD, &tmpv2).Reshape(n, 1); /* tmpv = (masks.GetSubmatrix(0, n - 1, i, i).GetHadamardProduct(x.GetSubmatrix(0, n - 1, j, j))) / sqn; */
                
                for(integer k = 0; k < n; k++)
                { /* ZY.SubmatrixAssignment(i * n, (i + 1) * n - 1, j, j, tmpv.Reshape(n1, n2).GetFFT2D(FFTW_FORWARD).Reshape(n, 1)); */
                    ZYptr[i * n + k + j * m].r = tmpv2ptr[k].r;
                    ZYptr[i * n + k + j * m].i = tmpv2ptr[k].i;
                }
            }
        }
        
        /*cD.SubmatrixAssignment(0, m - 1, 0, 0, ZY.GetTranspose().GetColNormsSquare().GetTranspose());
        cD = cD - b; */
        for(integer i = 0; i < m; i++)
        {
            for(integer j = 0; j < r; j++)
            {
                cDptr[i] += ZYptr[i + j * m].r * ZYptr[i + j * m].r + ZYptr[i + j * m].i * ZYptr[i + j * m].i;
            }
            cDptr[i] -= bptr[i];
        }
//        cD.SubmatrixAssignment(0, m - 1, 0, 0, ZY.GetTranspose().GetColNormsSquare().GetTranspose());
//        cD = cD - b;
        
        realdp result = cD.DotProduct(cD) / b.DotProduct(b);
        realdp tmpfn = x.DotProduct(x);
        result = (kappa == 0) ? result : result + kappa * tmpfn;
        x.AddToFields("cD", cD);
        x.AddToFields("ZY", ZY);
        return result;
        
        
//        Vector cD(m); cD.SetToZeros();
//        Vector ZY(m, r, "complex");
//        realdp sqn = sqrt(static_cast<realdp> (n1 * n2));
//        Vector tmpv;
//
//        for(integer i = 0; i < l; i++)
//        {
//            for(integer j = 0; j < r; j++)
//            {
//                tmpv = (masks.GetSubmatrix(0, n - 1, i, i).GetHadamardProduct(x.GetSubmatrix(0, n - 1, j, j))) / sqn;
//                ZY.SubmatrixAssignment(i * n, (i + 1) * n - 1, j, j, tmpv.Reshape(n1, n2).GetFFT2D(FFTW_FORWARD).Reshape(n, 1));
//            }
//        }
//        cD.SubmatrixAssignment(0, m - 1, 0, 0, ZY.GetTranspose().GetColNormsSquare().GetTranspose());
//        cD = cD - b;
//
//        realdp result = cD.DotProduct(cD) / b.DotProduct(b);
//        result = (kappa == 0) ? result : result + kappa * x.Fnorm() * x.Fnorm();
//        x.AddToFields("cD", cD);
//        x.AddToFields("ZY", ZY);
//        return result;
	};

	Vector &CSFRQPhaseRetrieval::EucGrad(const Variable &x, Vector *result) const
	{
        Vector cD = x.Field("cD"), ZY = x.Field("ZY");
        Vector tmpv(n, "complex"), tmpv2(n, "complex");
        Vector gfi(n, r, "complex");
        realdpcomplex *gfiptr = (realdpcomplex *) gfi.ObtainWriteEntireData();
        const realdpcomplex *masksptr = (realdpcomplex *) masks.ObtainReadData();
        realdpcomplex *tmpvptr = nullptr;
        result->SetToZeros();
        
        Vector DZY = cD.GetRealToComplex().GetDiagTimesM(ZY, GLOBAL::L);

        realdp sqn = sqrt(static_cast<realdp> (n1 * n2));
        for(integer i = 0; i < l; i++)
        {
            for(integer j = 0; j < r; j++)
            {
                /* tmpv = (DZY.GetSubmatrix(i * n, (i + 1) * n - 1, j, j).Reshape(n1, n2).GetFFT2D(FFTW_BACKWARD).Reshape(n, 1)) / sqn; */
                (DZY.GetSubmatrix(i * n, (i + 1) * n - 1, j, j).Reshape(n1, n2).FFT2D(FFTW_BACKWARD, &tmpv).Reshape(n, 1));
                tmpv.ScalarTimesThis(1.0 / sqn);
                tmpvptr = (realdpcomplex *) tmpv.ObtainWritePartialData();
                
                /* gfi.SubmatrixAssignment(0, n - 1, j, j, tmpv.GetHadamardProduct(masks.GetSubmatrix(0, n - 1, i, i).GetConj())); */
                for(integer k = 0; k < n; k++)
                {
                    gfiptr[k + j * n].r = tmpvptr[k].r * masksptr[k + i * n].r + tmpvptr[k].i * masksptr[k + i * n].i;
                    gfiptr[k + j * n].i = tmpvptr[k].i * masksptr[k + i * n].r - tmpvptr[k].r * masksptr[k + i * n].i;
                }
            }
            result->AlphaXaddThis(1, gfi);
        }
        result->ScalarTimesThis(4.0 / b.DotProduct(b));
        
        if(kappa != 0)
        {
            result->AlphaXaddThis(2.0 * kappa, x);
        }
        return *result;
	};

    Vector &CSFRQPhaseRetrieval::EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const
    {
        Vector cD = x.Field("cD"), ZY = x.Field("ZY");
        
        /*Zetax = Z * etax*/
        Vector Zetax(m, r, "complex");
        realdp sqn = sqrt(static_cast<realdp> (n1 * n2));
        Vector tmpv(n, "complex"), tmpv2(n, "complex");
        realdpcomplex *tmpvptr = (realdpcomplex *) tmpv.ObtainWriteEntireData();
        const realdpcomplex *masksptr = (realdpcomplex *) masks.ObtainReadData();
        const realdpcomplex *etaxptr = (realdpcomplex *) etax.ObtainReadData();
        
        for(integer i = 0; i < l; i++)
        {
            for(integer j = 0; j < r; j++)
            {
                /* tmpv = (masks.GetSubmatrix(0, n - 1, i, i).GetHadamardProduct(etax.GetSubmatrix(0, n - 1, j, j))) / sqn; */
                for(integer k = 0; k < n; k++)
                {
                    tmpvptr[k].r = (masksptr[k + i * n].r * etaxptr[k + j * n].r - masksptr[k + i * n].i * etaxptr[k + j * n].i) / sqn;
                    tmpvptr[k].i = (masksptr[k + i * n].i * etaxptr[k + j * n].r + masksptr[k + i * n].r * etaxptr[k + j * n].i) / sqn;
                }
                
                /* Zetax.SubmatrixAssignment(i * n, (i + 1) * n - 1, j, j, tmpv.Reshape(n1, n2).GetFFT2D(FFTW_FORWARD).Reshape(n, 1)); */
                tmpv.Reshape(n1, n2).FFT2D(FFTW_FORWARD, &tmpv2).Reshape(n, 1);
                Zetax.SubmatrixAssignment(i * n, (i + 1) * n - 1, j, j, tmpv2);
            }
        }
        
        Vector cdD = 2.0 * (Zetax.GetTranspose().GetColDotProducts(ZY.GetTranspose()));
        Vector dDZYaddDZeta = cdD.GetRealToComplex().GetDiagTimesM(ZY, GLOBAL::L) + cD.GetRealToComplex().GetDiagTimesM(Zetax, GLOBAL::L);
        Vector tmpi(n, r, "complex");
        realdpcomplex *tmpiptr = (realdpcomplex *) tmpi.ObtainWriteEntireData();
        result->SetToZeros();
        for(integer i = 0; i < l; i++)
        {
            for(integer j = 0; j < r; j++)
            {
                /*tmpv = (dDZYaddDZeta.GetSubmatrix(i * n, (i + 1) * n - 1, j, j).Reshape(n1, n2).GetFFT2D(FFTW_BACKWARD).Reshape(n, 1)) / sqn; */
                dDZYaddDZeta.GetSubmatrix(i * n, (i + 1) * n - 1, j, j).Reshape(n1, n2).FFT2D(FFTW_BACKWARD, &tmpv).Reshape(n, 1);
                tmpv.ScalarTimesThis(1.0 / sqn);
                
                /* tmpi.SubmatrixAssignment(0, n - 1, j, j, tmpv.GetHadamardProduct(masks.GetSubmatrix(0, n - 1, i, i).GetConj())); */
                for(integer k = 0; k < n; k++)
                {
                    tmpiptr[k + j * n].r = tmpvptr[k].r * masksptr[k + i * n].r + tmpvptr[k].i * masksptr[k + i * n].i;
                    tmpiptr[k + j * n].i = tmpvptr[k].i * masksptr[k + i * n].r - tmpvptr[k].r * masksptr[k + i * n].i;
                }
            }
            result->AlphaXaddThis(1, tmpi);
        }
        result->ScalarTimesThis(4.0 / b.DotProduct(b));
        result->AlphaXaddThis(2.0 * kappa, etax); /*result * (4.0 / b.DotProduct(b)) + 2.0 * kappa * etax; */
        return *result;
    };

}; /*end of ROPTLIB namespace*/
#endif
