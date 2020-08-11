
#include "Manifolds/FixedRankQ2F.h"

/*Define the namespace*/
namespace ROPTLIB{

	FixedRankQ2F::FixedRankQ2F(integer inm, integer inn, integer inr) : ProductManifold(2,
		new Euclidean(inm, inr), static_cast<integer> (1), new Euclidean(inn, inr), static_cast<integer> (1))
	{
		m = inm;
		n = inn;
		r = inr;
		name.assign("Fixed-rank manifold by 2-factor form");
		IsIntrApproach = true;
        Vector F1(m, r), F2(n, r);
        Vector Prod(2, &F1, 1, &F2, 1);
        EMPTYEXTR = Prod;
		EMPTYINTR = Vector ((m + n - r) * r);
	};

    void FixedRankQ2F::CheckParams(void) const
    {
        Manifold::CheckParams();
        printf("%s PARAMETERS:\n", name.c_str());
        printf("row (m)       :%15d,\t", m);
        printf("column (n)    :%15d,\t", n);
        printf("rank (r)      :%15d\n", r);
    };

    Variable FixedRankQ2F::RandominManifold(void) const
    {
        Variable result(EMPTYEXTR); result.RandGaussian();
        return result;
    };

	FixedRankQ2F::~FixedRankQ2F()
	{
		for (integer i = 0; i < numoftypes; i++)
		{
			delete manifolds[i];
		}
	};

    Vector &FixedRankQ2F::Projection(const Variable &x, const Vector &etax, Vector *result) const
    {
        if(IsIntrApproach)
        {
            *result = etax;
            return *result;
        }
        
        return ExtrProjection(x, etax, result);
    };

    Vector &FixedRankQ2F::ExtrProjection(const Variable &x, const Vector &etax, Vector *result) const
    {
        GenerateFieldsX(x);
        Vector LG = x.GetElement(0).Field("_L"), LH = x.GetElement(1).Field("_L");
        
        /*tmp1 = (H^T H)^{-1} (H^T Z)*/
        Vector inetax = etax;
        Vector tmp1(r, r); tmp1.AlphaABaddBetaThis(1, x.GetElement(1), GLOBAL::T, inetax.GetElement(1), GLOBAL::N, 0);
        tmp1 = tmp1.TriangleLinSol(LH, GLOBAL::N);
        tmp1 = tmp1.TriangleLinSol(LH, GLOBAL::T);
        
        /*tmp2 = (G^T G)^{-1} (G^T W)*/
        Vector tmp2(r, r); tmp2.AlphaABaddBetaThis(1, x.GetElement(0), GLOBAL::T, inetax.GetElement(0), GLOBAL::N, 0);
        tmp2 = tmp2.TriangleLinSol(LG, GLOBAL::N);
        tmp2 = tmp2.TriangleLinSol(LG, GLOBAL::T);

        Vector Lambda = (tmp1.GetTranspose() - tmp2) / 2;
        result->NewMemoryOnWrite();
        result->GetElement(0) = inetax.GetElement(0);
        result->GetElement(0).AlphaABaddBetaThis(1, x.GetElement(0), GLOBAL::N, Lambda, GLOBAL::N, 1); /* result.GetElement(0) = W + G * Lambda; */
        result->GetElement(1) = inetax.GetElement(1);
        result->GetElement(1).AlphaABaddBetaThis(-1, x.GetElement(1), GLOBAL::N, Lambda, GLOBAL::T, 1); /* result.GetElement(1) = Z - H * Lambda.GetTranspose(); */
        return *result;
    };

    realdp FixedRankQ2F::Metric(const Variable &x, const Vector &etax, const Vector &xix) const
    {
        if(IsIntrApproach)
            return Manifold::Metric(x, etax, xix);
        
        Vector G = x.GetElement(0), H = x.GetElement(1);
        Vector etaxG = etax.GetElement(0), etaxH = etax.GetElement(1);
        Vector xixG = xix.GetElement(0), xixH = xix.GetElement(1);
        
        realdp result = 0;
        Vector tmp1(r, r); tmp1.AlphaABaddBetaThis(1, G, GLOBAL::T, G, GLOBAL::N, 0);
        Vector tmp2(r, r); tmp2.AlphaABaddBetaThis(1, xixH, GLOBAL::T, etaxH, GLOBAL::N, 0);
        result += tmp1.DotProduct(tmp2);
        tmp1.AlphaABaddBetaThis(1, H, GLOBAL::T, H, GLOBAL::N, 0); tmp2.AlphaABaddBetaThis(1, etaxG, GLOBAL::T, xixG, GLOBAL::N, 0);
        result += tmp1.DotProduct(tmp2);
        tmp1.AlphaABaddBetaThis(1, H, GLOBAL::T, etaxH, GLOBAL::N, 0); tmp2.AlphaABaddBetaThis(1, xixG, GLOBAL::T, G, GLOBAL::N, 0);
        result += tmp1.DotProduct(tmp2);
        tmp1.AlphaABaddBetaThis(1, xixH, GLOBAL::T, H, GLOBAL::N, 0); tmp2.AlphaABaddBetaThis(1, G, GLOBAL::T, etaxG, GLOBAL::N, 0);
        result += tmp1.DotProduct(tmp2);
        
        return result;
    };

	Vector &FixedRankQ2F::ObtainIntr(const Variable &x, const Vector &etax, Vector *result) const
	{
        Variable G = x.GetElement(0);
        Variable H = x.GetElement(1);
        if(!G.FieldsExist("_HHR"))
        {
            G.HHRDecom();
        }
        if(!H.FieldsExist("_HHR"))
        {
            H.HHRDecom();
        }
        Vector etaxG = etax.GetElement(0);
        Vector etaxH = etax.GetElement(1);

        /*for G component*/
        Vector Gtmp = etaxG.HHRMtp(G.Field("_HHR"), G.Field("_tau"), GLOBAL::T, GLOBAL::L);
        Vector GHHR = G.Field("_HHR");
        const realdp *GHHRptr = GHHR.ObtainReadData();
        realdp *Gtmpptr = Gtmp.ObtainWritePartialData();
        for (integer i = 0; i < r; i++)
        {
            if(GHHRptr[i + m * i] < 0)
                scal_(&r, &GLOBAL::DNONE, Gtmpptr + i, &m);
        }
        
        if(!G.FieldsExist("_LG"))
        {
            Vector LG(r, r);
            realdp *LGptr = LG.ObtainWriteEntireData();
            realdp sign = 0;
            for (integer i = 0; i < r; i++)
            {
                for (integer j = 0; j < i; j++)
                {
                    LGptr[j + i * r] = 0;
                }
                sign = ((GHHRptr[i + m * i] >= 0) ? 1 : -1);
                for (integer j = i; j < r; j++)
                {
                    LGptr[j + i * r] = GHHRptr[i + m * j] * sign;
                }
            }
            G.AddToFields("_LG", LG);
        }
        Vector LG = G.Field("_LG");
        
        /*for H component*/
        Vector Htmp = etaxH.HHRMtp(H.Field("_HHR"), H.Field("_tau"), GLOBAL::T, GLOBAL::L);
        Vector HHHR = H.Field("_HHR");
        realdp *HHHRptr = HHHR.ObtainWritePartialData();
        realdp *Htmpptr = Htmp.ObtainWritePartialData();
        for (integer i = 0; i < r; i++)
        {
            if(HHHRptr[i + n * i] < 0)
                scal_(&r, &GLOBAL::DNONE, Htmpptr + i, &n);
        }
        if(!H.FieldsExist("_LH"))
        {
            Vector LH(r, r);
            realdp *LHptr = LH.ObtainWriteEntireData();
            realdp sign = 0;
            for (integer i = 0; i < r; i++)
            {
                for (integer j = 0; j < i; j++)
                {
                    LHptr[j + i * r] = 0;
                }
                sign = ((HHHRptr[i + n * i] >= 0) ? 1 : -1);
                for (integer j = i; j < r; j++)
                {
                    LHptr[j + i * r] = HHHRptr[i + n * j] * sign;
                }
            }
            H.AddToFields("_LH", LH);
        }
        Vector LH = H.Field("_LH");

        Vector GtmpLH = Gtmp * LH, HtmpLG = Htmp * LG;
        const realdp *GtmpLHptr = GtmpLH.ObtainReadData();
        const realdp *HtmpLGptr = HtmpLG.ObtainReadData();
        
        realdp *resultptr = result->ObtainWriteEntireData();
        integer idx = 0;
        for(integer i = 0; i < r; i++)
        {
            for(integer j = 0; j < r; j++)
            {
                resultptr[idx] = (GtmpLHptr[j + i * m] + HtmpLGptr[i + j * n]);
                idx++;
            }
        }
        
        for (integer i = 0; i < r; i++)
        {
            for (integer j = r; j < m; j++)
            {
                resultptr[idx] = GtmpLHptr[j + i * m];
                idx++;
            }
        }
        for (integer i = 0; i < r; i++)
        {
            for (integer j = r; j < n; j++)
            {
                resultptr[idx] = HtmpLGptr[j + i * n];
                idx++;
            }
        }
        
        return *result;
	};

	Vector &FixedRankQ2F::ObtainExtr(const Variable &x, const Vector &intretax, Vector *result) const
	{
        Variable G = x.GetElement(0);
        Variable H = x.GetElement(1);
        if(!G.FieldsExist("_HHR"))
        {
            G.HHRDecom();
        }
        if(!H.FieldsExist("_HHR"))
        {
            H.HHRDecom();
        }
        
        const realdp *intretaxptr = intretax.ObtainReadData();
        
        result->NewMemoryOnWrite();
        Vector etaxG = result->GetElement(0);
        Vector etaxH = result->GetElement(1);
        realdp *etaxGptr = etaxG.ObtainWriteEntireData();
        realdp *etaxHptr = etaxH.ObtainWriteEntireData();

        integer idx = 0;
        for (integer i = 0; i < r; i++)
        {
            for (integer j = 0; j < r; j++)
            {
                etaxGptr[j + i * m] = intretaxptr[idx] / 2;
                etaxHptr[i + j * n] = intretaxptr[idx] / 2;
                idx++;
            }
        }
        for (integer i = 0; i < r; i++)
        {
            for (integer j = r; j < m; j++)
            {
                etaxGptr[j + i * m] = intretaxptr[idx];
                idx++;
            }
        }
        for (integer i = 0; i < r; i++)
        {
            for (integer j = r; j < n; j++)
            {
                etaxHptr[j + i * n] = intretaxptr[idx];
                idx++;
            }
        }

        /*for etaG*/
        Vector GHHR = G.Field("_HHR");
        realdp *GHHRptr = GHHR.ObtainWritePartialData();
        for (integer i = 0; i < r; i++)
        {
            if(GHHRptr[i + m * i] < 0)
                scal_(&r, &GLOBAL::DNONE, etaxGptr + i, &m);
        }
        etaxG = etaxG.HHRMtp(G.Field("_HHR"), G.Field("_tau"), GLOBAL::N, GLOBAL::L);
        
        Vector LG(r, r);
        realdp *LGptr = LG.ObtainWriteEntireData();
        realdp sign = 0;

        for (integer i = 0; i < r; i++)
        {
            for (integer j = 0; j < i; j++)
            {
                LGptr[j + i * r] = 0;
            }
            sign = ((GHHRptr[i + m * i] >= 0) ? 1 : -1);
            for (integer j = i; j < r; j++)
            {
                LGptr[j + i * r] = GHHRptr[i + m * j] * sign;
            }
        }

        /*for etaH*/
        Vector HHHR = H.Field("_HHR");
        realdp *HHHRptr = HHHR.ObtainWritePartialData();
        for (integer i = 0; i < r; i++)
        {
            if(HHHRptr[i + n * i] < 0)
                scal_(&r, &GLOBAL::DNONE, etaxHptr + i, &n);
        }
        etaxH = etaxH.HHRMtp(H.Field("_HHR"), H.Field("_tau"), GLOBAL::N, GLOBAL::L);
        
        Vector LH(r, r);
        realdp *LHptr = LH.ObtainWriteEntireData();

        for (integer i = 0; i < r; i++)
        {
            for (integer j = 0; j < i; j++)
            {
                LHptr[j + i * r] = 0;
            }
            sign = ((HHHRptr[i + n * i] >= 0) ? 1 : -1);
            for (integer j = i; j < r; j++)
            {
                LHptr[j + i * r] = HHHRptr[i + n * j] * sign;
            }
        }

        result->GetElement(0) = etaxG.GetTranspose().TriangleLinSol(LH, GLOBAL::T).GetTranspose();
        
        result->GetElement(1) = etaxH.GetTranspose().TriangleLinSol(LG, GLOBAL::T).GetTranspose();

        return *result;
	};

    Vector &FixedRankQ2F::coTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
    {
        printf("warning:FixedRankQ2F::coTangentVector has not been done!\n");
        
        return MultiManifolds::coTangentVector(x, etax, y, xiy, result);
    };

    Vector &FixedRankQ2F::VectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const
    {
        if (HasHHR)
            return LCVectorTransport(x, etax, y, xix, result);
        
        if(IsIntrApproach)
        {
            *result = xix;
            return *result;
        }
        Vector inxix(EMPTYINTR); ObtainIntr(x, xix, &inxix);
        return ObtainExtr(y, inxix, result);
    };

    Vector &FixedRankQ2F::InverseVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
    {
        if (HasHHR)
            return LCInverseVectorTransport(x, etax, y, xiy, result);
        
        if(IsIntrApproach)
        {
            *result = xiy;
            return *result;
        }
        
        Vector inxiy(EMPTYINTR); ObtainIntr(x, xiy, &inxiy);
        return ObtainExtr(x, inxiy, result);
    };

    LinearOPE &FixedRankQ2F::TranHInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, LinearOPE *result) const
    {
        if (HasHHR)
            return LCTranHInvTran(x, etax, y, Hx, result);
        
        if(IsIntrApproach)
        {
            *result = Hx;
            return *result;
        }
        
        return MultiManifolds::TranHInvTran(x, etax, y, Hx, result);
    };

    LinearOPE &FixedRankQ2F::HaddScaledRank1OPE(const Variable &x, const LinearOPE &Hx, realdp scalar, const Vector &etax, const Vector &xix, LinearOPE *result) const
    {
        if(IsIntrApproach)
            return Manifold::HaddScaledRank1OPE(x, Hx, scalar, etax, xix, result);
        
        return MultiManifolds::HaddScaledRank1OPE(x, Hx, scalar, etax, xix, result);
    };

	Variable &FixedRankQ2F::Retraction(const Variable &x, const Vector &etax, Variable *result) const
	{
        Vector exetax(EMPTYEXTR);
        if(IsIntrApproach)
            ObtainExtr(x, etax, &exetax);
        else
            exetax = etax;

        /*Simple retraction*/
        ProductManifold::Retraction(x, exetax, result);

        GenerateFieldsX(*result);

        return *result;
	};

    Variable FixedRankQ2F::ProjRetraction(const Variable &x, const Vector &etax, Variable *result) const
    { /*GH^T + dG H^T + G dH^T = [G dG] [I, I; I, 0] * [H dH]^T*/
        Vector exetax(EMPTYEXTR);
        if(IsIntrApproach)
            ObtainExtr(x, etax, &exetax);
        else
            exetax = etax;
        
        /*Retraction by the projection in its embedding space. It is by svd*/
        Vector G = x.GetElement(0), H = x.GetElement(1);
        Vector dG = exetax.GetElement(0), dH = exetax.GetElement(1);
        Vector tmp1(m, 2 * r), tmp2(n, 2 * r);
        realdp *tmp1ptr = tmp1.ObtainWriteEntireData();
        realdp *tmp2ptr = tmp2.ObtainWriteEntireData();
        const realdp *Gptr = G.ObtainReadData(), *Hptr = H.ObtainReadData();
        const realdp *dGptr = dG.ObtainReadData(), *dHptr = dH.ObtainReadData();
        integer length = m * r;
        copy_(&length, const_cast<realdp *> (Gptr), &GLOBAL::IONE, tmp1ptr, &GLOBAL::IONE);
        copy_(&length, const_cast<realdp *> (dGptr), &GLOBAL::IONE, tmp1ptr + length, &GLOBAL::IONE); /*tmp1 = [G dG]*/
        length = n * r;
        copy_(&length, const_cast<realdp *> (Hptr), &GLOBAL::IONE, tmp2ptr, &GLOBAL::IONE);
        copy_(&length, const_cast<realdp *> (dHptr), &GLOBAL::IONE, tmp2ptr + length, &GLOBAL::IONE); /*tmp2 = [H dH]*/
        
        integer *ir = new integer[6 * r];
        integer *jc = ir + 3 * r;
        realdp *vals = new realdp [3 * r];
        for(integer i = 0; i < r; i++)
        {
            ir[i] = i;      ir[i + r] = i;      ir[i + 2 * r] = i + r;
            jc[i] = i;      jc[i + r] = i + r;  jc[i + 2 * r] = i;
            vals[i] = 1;    vals[i + r] = 1;    vals[i + 2 * r] = 1;
        }
        
        SparseMatrix SM(2 * r, 2 * r, ir, jc, vals, 3 * r); /* SM = [I, I; I, 0] */
        
        tmp1.QRDecom(); tmp2.QRDecom();
        Vector tmp1Q = tmp1.Field("_Q"), tmp1R = tmp1.Field("_R"); /*[G dG] = tmp1Q * tmp1R is a QR decomposition*/
        Vector tmp2Q = tmp2.Field("_Q"), tmp2R = tmp2.Field("_R"); /*[H dH] = tmp2Q * tmp2R is a QR decomposition*/
        
        Vector M = tmp1R * SM * tmp2R.GetTranspose();
        M.SVDDecom();
        Vector U = M.Field("_U"), S = M.Field("_S"), Vt = M.Field("_Vt");

        realdp *Sptr = S.ObtainWritePartialData();
        for(integer i = 0; i < 2 * r; i++)
        {
            Sptr[i] = std::sqrt(Sptr[i]);
        }
        
        S = S.GetSubmatrix(0, r - 1, 0, 0);
        
        result->NewMemoryOnWrite();
        result->GetElement(0) = S.GetDiagTimesM(tmp1Q * U.GetSubmatrix(0, 2 * r - 1, 0, r - 1), GLOBAL::R);
        result->GetElement(1) = S.GetDiagTimesM(tmp2Q * Vt.GetTranspose().GetSubmatrix(0, 2 * r - 1, 0, r - 1), GLOBAL::R);
        
        GenerateFieldsX(*result);

        delete [] ir;
        delete [] vals;
        return *result;
    };

	Vector &FixedRankQ2F::DiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir) const
	{
        realdp nxix = std::sqrt(Metric(x, xix, xix));
        
        Vector exxix(EMPTYEXTR);
        
        if(IsIntrApproach)
            ObtainExtr(x, xix, &exxix);
        else
            exxix = xix;
        
        Vector exresult(EMPTYEXTR);
        ExtrProjection(y, exxix, &exresult);
        
        if(IsIntrApproach)
            ObtainIntr(y, exresult, result);
        else
            *result = exresult;
        
        if (IsEtaXiSameDir && HasHHR)
        {
            Vector beta(3);
            realdp *betaptr = beta.ObtainWriteEntireData();
            realdp EtatoXi = std::sqrt(Metric(x, etax, etax)) / nxix;
            betaptr[0] = std::sqrt(Metric(x, etax, etax) / Metric(y, *result, *result)) / EtatoXi;
            betaptr[1] = Metric(x, etax, etax);
            betaptr[2] = Metric(y, *result, *result) * EtatoXi * EtatoXi;
            etax.AddToFields("beta", beta);
            
            if (HasHHR)
            {
                etax.AddToFields("betaTReta", (*result) * (betaptr[0] * EtatoXi));
            }
        }
        
        return *result;
	};

	Vector &FixedRankQ2F::EucGradToGrad(const Variable &x, const Vector &egf, const Problem *prob, Vector *result)  const
	{
        if(prob->GetUseHess())
        {
            /*The copy on write is necessary. The reason is that the egf may be from a component in a product of elements.
            Therefore, if CopyOnWrite is not used, then the attached data in x and the product of elements share the same
            memory. This may cause an issue: if the product of elements are released before the attached data in x, then
            release the attached data in x would attempt to delete memory that has been released. This is an error!*/
            Vector Sharedegf(egf);
            Sharedegf.CopyOnWrite();
            x.AddToFields("EGrad", Sharedegf);
        }
        
        result->NewMemoryOnWrite();
        
        Vector G = x.GetElement(0), H = x.GetElement(1);
        if(!G.FieldsExist("_LG"))
        {
            Vector tmp = G.GetTranspose() * G;
            tmp.CholDecom();
            x.GetElement(0).AddToFields("_L", tmp.Field("_L"));
        }
        if(!H.FieldsExist("_LH"))
        {
            Vector tmp = H.GetTranspose() * H;
            tmp.CholDecom();
            x.GetElement(1).AddToFields("_L", tmp.Field("_L"));
        }
        
        Vector egfG = egf.GetElement(0), egfH = egf.GetElement(1);
        Vector LG = x.GetElement(0).Field("_L"), LH = x.GetElement(1).Field("_L");
        
        Vector tmpG = egfG.GetTranspose().TriangleLinSol(LH, GLOBAL::N).TriangleLinSol(LH, GLOBAL::T).GetTranspose();
        result->GetElement(0) = tmpG - 0.5 * G * (G.GetTranspose() * tmpG).TriangleLinSol(LG, GLOBAL::N).TriangleLinSol(LG, GLOBAL::T);
        
        Vector tmpH = egfH.GetTranspose().TriangleLinSol(LG, GLOBAL::N).TriangleLinSol(LG, GLOBAL::T).GetTranspose();
        result->GetElement(1) = tmpH - 0.5 * H * (H.GetTranspose() * tmpH).TriangleLinSol(LH, GLOBAL::N).TriangleLinSol(LH, GLOBAL::T);
        
        return *result;
	};

	Vector &FixedRankQ2F::EucHvToHv(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const
	{
        Vector G = x.GetElement(0), H = x.GetElement(1);
        Vector dG = etax.GetElement(0), dH = etax.GetElement(1);
        Vector EGrad = x.Field("EGrad");
        Vector EG_dG = EGrad.GetElement(0), EG_dH = EGrad.GetElement(1);
        Vector EH_dG = exix.GetElement(0), EH_dH = exix.GetElement(1);
        
        result->NewMemoryOnWrite();
        
        Vector LG = x.GetElement(0).Field("_L"), LH = x.GetElement(1).Field("_L");
        
        Vector tmpG = EH_dG.GetTranspose().TriangleLinSol(LH, GLOBAL::N).TriangleLinSol(LH, GLOBAL::T).GetTranspose();
        tmpG = tmpG - 0.5 * G * (G.GetTranspose() * tmpG).TriangleLinSol(LG, GLOBAL::N).TriangleLinSol(LG, GLOBAL::T);
        
        Vector tmpH = EH_dH.GetTranspose().TriangleLinSol(LG, GLOBAL::N).TriangleLinSol(LG, GLOBAL::T).GetTranspose();
        tmpH = tmpH - 0.5 * H * (H.GetTranspose() * tmpH).TriangleLinSol(LH, GLOBAL::N).TriangleLinSol(LH, GLOBAL::T);
        
        /*tmp1 = EG_dG * (H^T H)^{-1} * (H^T dH) * (H^T H)^{-1} */
        Vector tmp1 = EG_dG * ((dH.GetTranspose() * H).TriangleLinSol(LH, GLOBAL::N).TriangleLinSol(LH, GLOBAL::T).GetTranspose().TriangleLinSol(LH, GLOBAL::N).TriangleLinSol(LH, GLOBAL::T));
        
        /*tmp1 = (G (G^T G)^{-1} G^T - I) * EG_dG * (H^T H)^{-1} * (H^T dH) * (H^T H)^{-1} */
        tmp1 = G * (G.GetTranspose() * tmp1).TriangleLinSol(LG, GLOBAL::N).TriangleLinSol(LG, GLOBAL::T) - tmp1;
        
        /* tmp1 = (G (G^T G)^{-1} G^T - I) * EG_dG * (H^T H)^{-1} * (H^T dH) * (H^T H)^{-1} - 0.5 * G (G^T G)^{-1} (EG_dH^T dH) * (H^T H)^{-1}*/
        tmp1 = tmp1 - 0.5 * G * ((dH.GetTranspose() * EG_dH).TriangleLinSol(LH, GLOBAL::N).TriangleLinSol(LH, GLOBAL::T).GetTranspose().TriangleLinSol(LG, GLOBAL::N).TriangleLinSol(LG, GLOBAL::T));
        result->GetElement(0) = tmpG + tmp1;
        
        /*tmp2 = EG_dH * (G^T G)^{-1} * (G^T dG) * (G^T G)^{-1} */
        Vector tmp2 = EG_dH * ((dG.GetTranspose() * G).TriangleLinSol(LG, GLOBAL::N).TriangleLinSol(LG, GLOBAL::T).GetTranspose().TriangleLinSol(LG, GLOBAL::N).TriangleLinSol(LG, GLOBAL::T));
        
        /*tmp2 = (H (H^T H)^{-1} H^T - I) * EG_dH * (G^T G)^{-1} * (G^T dG) * (G^T G)^{-1} */
        tmp2 = H * (H.GetTranspose() * tmp2).TriangleLinSol(LH, GLOBAL::N).TriangleLinSol(LH, GLOBAL::T) - tmp2;
        
        /* tmp2 = (H (H^T H)^{-1} H^T - I) * EG_dH * (G^T G)^{-1} * (G^T dG) * (G^T G)^{-1} - 0.5 * H (H^T H)^{-1} (EG_dG^T dG) * (G^T G)^{-1}*/
        tmp2 = tmp2 - 0.5 * H * ((dG.GetTranspose() * EG_dG).TriangleLinSol(LG, GLOBAL::N).TriangleLinSol(LG, GLOBAL::T).GetTranspose().TriangleLinSol(LH, GLOBAL::N).TriangleLinSol(LH, GLOBAL::T));
        result->GetElement(1) = tmpH + tmp2;
        
        return ExtrProjection(x, *result, result);
	};

    void FixedRankQ2F::GenerateFieldsX(const Variable &x) const
    {
        Variable G = x.GetElement(0);
        Variable H = x.GetElement(1);
        if(!G.FieldsExist("_LG"))
        {
            Vector tmp = G.GetTranspose() * G;
            tmp.CholDecom();
            x.GetElement(0).AddToFields("_L", tmp.Field("_L"));
        }
        
        if(!H.FieldsExist("_LH"))
        {
            Vector tmp = H.GetTranspose() * H;
            tmp.CholDecom();
            x.GetElement(1).AddToFields("_L", tmp.Field("_L"));
        }
    };
}; /*end of ROPTLIB namespace*/
