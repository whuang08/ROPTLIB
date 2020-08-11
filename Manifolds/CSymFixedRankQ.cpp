
#include "Manifolds/CSymFixedRankQ.h"

/*Define the namespace*/
namespace ROPTLIB{

	CSymFixedRankQ::CSymFixedRankQ(integer r, integer c)
	{
		n = r;
		p = c;
		IsIntrApproach = true;
		HasHHR = false;
		name.assign("Hermitian fixed rank manifold by quotient representation");
		IntrinsicDim = 2 * r * c - c * c;
		ExtrinsicDim = 2 * r * c;
        metric = CSFRANKQEUC;
        VecTran = CSFRANKQPARALLELIZATION;
        EMPTYEXTR = Vector(r, c, "complex");
        EMPTYINTR = Vector(IntrinsicDim);
	};

	CSymFixedRankQ::~CSymFixedRankQ(void)
	{
	};

    void CSymFixedRankQ::SetParams(PARAMSMAP params)
    {
        Manifold::SetParams(params);
        PARAMSMAP::iterator iter;
        for (iter = params.begin(); iter != params.end(); iter++)
        {
            if (iter->first == static_cast<std::string> ("ParamSet"))
            {
                switch (static_cast<integer> (iter->second))
                {
                case 1:
                    ChooseParamsSet1();
                    break;
                case 2:
                    ChooseParamsSet2();
                    break;
                case 3:
                    ChooseParamsSet3();
                    break;
                case 4:
                    ChooseParamsSet4();
                    break;
                default:
                    printf("warning: invalid index for parameter set! use the default one.");
                    ChooseParamsSet1();
                    break;
                }
            }
        }
    };

    void CSymFixedRankQ::ChooseParamsSet1(void)
    {
        metric = CSFRANKQEUC;
        VecTran = CSFRANKQPARALLELIZATION;
    };

    void CSymFixedRankQ::ChooseParamsSet2(void)
    {
        metric = CSFRANKQHGZ;
        VecTran = CSFRANKQPARALLELIZATION;
    };

    void CSymFixedRankQ::ChooseParamsSet3(void)
    {
        metric = CSFRANKQEUC;
        VecTran = CSFRANKQPROJECTION;
    };

    void CSymFixedRankQ::ChooseParamsSet4(void)
    {
        metric = CSFRANKQHGZ;
        VecTran = CSFRANKQPROJECTION;
    };

    Variable CSymFixedRankQ::RandominManifold(void) const
    {
        Variable result(EMPTYEXTR); result.RandGaussian();
        return result;
    };

    void CSymFixedRankQ::CheckParams(void) const
    {
        std::string SFRANKQMetricnames[CSFRANKQMETRICLENGTH] = { "CSFRANKQEUC", "SFRANKQHGZ" };
        std::string SFRANKQVTnames[CSFRANKQVECTORTRANSPORTLENGTH] = { "CSFRANKQPARALLELIZATION", "CSFRANKQPROJECTION" };
        Manifold::CheckParams();
        printf("%s PARAMETERS:\n", name.c_str());
        if (p == 1)
        {
            printf("n           :%15d,\t", n);
            printf("metric      :%15s\n", SFRANKQMetricnames[metric].c_str());
            printf("VecTran     :%15s\n", SFRANKQVTnames[VecTran].c_str());
        }
        else
        {
            printf("n           :%15d,\t", n);
            printf("p           :%15d\n", p);
            printf("metric      :%15s\t", SFRANKQMetricnames[metric].c_str());
            printf("VecTran     :%15s\n", SFRANKQVTnames[VecTran].c_str());
        }
    };

    realdp CSymFixedRankQ::Metric(const Variable &x, const Vector &etax, const Vector &xix) const
    {
        if (IsIntrApproach)
            return Manifold::Metric(x, etax, xix);
        
        if(metric == CSFRANKQEUC)
        {/*2 \trace(x^H etax x^H xix + x^H x etax^H xix)*/
            
            Vector tmp1(p, p, "complex"); tmp1.AlphaABaddBetaThis(GLOBAL::ZONE, etax, GLOBAL::C, x, GLOBAL::N, GLOBAL::ZZERO);
            Vector tmp2(p, p, "complex"); tmp2.AlphaABaddBetaThis(GLOBAL::ZONE, x, GLOBAL::C, xix, GLOBAL::N, GLOBAL::ZZERO);
            realdp result = tmp1.DotProduct(tmp2);
            tmp1.AlphaABaddBetaThis(GLOBAL::ZONE, x, GLOBAL::C, x, GLOBAL::N, GLOBAL::ZZERO);
            tmp2.AlphaABaddBetaThis(GLOBAL::ZONE, etax, GLOBAL::C, xix, GLOBAL::N, GLOBAL::ZZERO);
            result += tmp1.DotProduct(tmp2);
            return 2.0 * result; /* 2.0 * ((etax.GetTranspose() * x).DotProduct(x.GetTranspose() * xix) + (x.GetTranspose() * x).DotProduct(etax.GetTranspose() * xix)); */
        }
        
        if(metric == CSFRANKQHGZ)
        { /*\trace(x^H x etax^H xix)*/
            Vector tmp1(p, p, "complex"); tmp1.AlphaABaddBetaThis(GLOBAL::ZONE, x, GLOBAL::C, x, GLOBAL::N, GLOBAL::ZZERO);
            Vector tmp2(p, p, "complex"); tmp2.AlphaABaddBetaThis(GLOBAL::ZONE, etax, GLOBAL::C, xix, GLOBAL::N, GLOBAL::ZZERO);
            return tmp1.DotProduct(tmp2); /*(x.GetTranspose() * x).DotProduct(etax.GetTranspose() * xix); */
        }
        
        printf("Warning: this metric has not been done!\n");
        return 0;
    };

    Vector &CSymFixedRankQ::Projection(const Variable &x, const Vector &etax, Vector *result) const
    {
        if (IsIntrApproach)
            return IntrProjection(x, etax, result);
        
        return ExtrProjection(x, etax, result);
    };

    Vector &CSymFixedRankQ::ExtrProjection(const Variable &x, const Vector &etax, Vector *result) const
    {
        Vector XTX(p, p, "complex"); XTX.AlphaABaddBetaThis(GLOBAL::ZONE, x, GLOBAL::C, x, GLOBAL::N, GLOBAL::ZZERO);
        Vector XTV(p, p, "complex"); XTV.AlphaABaddBetaThis(GLOBAL::ZONE, x, GLOBAL::C, etax, GLOBAL::N, GLOBAL::ZZERO);
        XTV = XTX % XTV;
        
        *result = etax;
        XTX = (XTV - XTV.GetTranspose()) / 2.0;
        result->AlphaABaddBetaThis(GLOBAL::ZNONE, x, GLOBAL::N, XTX, GLOBAL::N, GLOBAL::ZONE); /* result = etax - x * ((XTV - XTV.GetTranspose()) / 2.0);*/
        return *result;
    };

    Vector &CSymFixedRankQ::ObtainIntr(const Variable &x, const Vector &etax, Vector *result) const
    {
        if (metric == CSFRANKQEUC)
            return ObtainIntrEUC(x, etax, result);
        
        if (metric == CSFRANKQHGZ)
            return ObtainIntrHGZ(x, etax, result);
            
        printf("warning: CSymFixedRankQ::ObtainIntr has not been done for this metric!\n");
        return Manifold::ObtainIntr(x, etax, result);
    };

    Vector &CSymFixedRankQ::ObtainExtr(const Variable &x, const Vector &intretax, Vector *result) const
    {
        if (metric == CSFRANKQEUC)
            return ObtainExtrEUC(x, intretax, result);
        
        if (metric == CSFRANKQHGZ)
            return ObtainExtrHGZ(x, intretax, result);
        
        printf("warning: CSymFixedRankQ::ObtainExtr has not been done for this metric!\n");
        return Manifold::ObtainExtr(x, intretax, result);
    };

    Vector &CSymFixedRankQ::coTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
    {
        printf("warning:CSymFixedRankQ::coTangentVector has not been done!\n");
        return Manifold::coTangentVector(x, etax, y, xiy, result);
    };

    Vector &CSymFixedRankQ::VectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const
    {
        if (VecTran == CSFRANKQPARALLELIZATION && !HasHHR)
        {
            return Manifold::VectorTransport(x, etax, y, xix, result);
        }

        if (VecTran == CSFRANKQPROJECTION && !HasHHR)
        {
            return VectorTransportProj(x, etax, y, xix, result);
        }

        if (HasHHR)
            return LCVectorTransport(x, etax, y, xix, result);

        printf("Error: CSymFixedRankQ::VectorTransport has not been done!\n");
        return Manifold::VectorTransport(x, etax, y, xix, result);
    };

    Vector &CSymFixedRankQ::InverseVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result)  const
    {
        if (VecTran == CSFRANKQPARALLELIZATION && !HasHHR)
        {
            return Manifold::InverseVectorTransport(x, etax, y, xiy, result);
        }

        if (VecTran == CSFRANKQPROJECTION && !HasHHR)
        {
            printf("Warning: Inverse Vector transport by projection has not been done!\n");
            return Manifold::InverseVectorTransport(x, etax, y, xiy, result);
        }

        if (HasHHR)
            return LCInverseVectorTransport(x, etax, y, xiy, result);

        printf("Error: CSymFixedRankQ::InverseVectorTransport has not been done!\n");
        return Manifold::InverseVectorTransport(x, etax, y, xiy, result);
    };

    LinearOPE &CSymFixedRankQ::TranHInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, LinearOPE *result) const
    {
        if (VecTran == CSFRANKQPARALLELIZATION && !HasHHR)
        {
            return Manifold::TranHInvTran(x, etax, y, Hx, result);
        }

        if (VecTran == CSFRANKQPROJECTION && !HasHHR)
        {
            printf("Warning: CSymFixedRankQ::TranHInvTran: Transport a linear operator using vector transport by projection has not been done!\n");
            return Manifold::TranHInvTran(x, etax, y, Hx, result);
        }

        if (HasHHR)
            return LCVectorTransport(x, etax, y, Hx, result);

        printf("Error: CSymFixedRankQ::TranHInvTran has not been done!\n");
        return Manifold::TranHInvTran(x, etax, y, Hx, result);
    };

    Variable &CSymFixedRankQ::Retraction(const Variable &x, const Vector &etax, Variable *result) const
    {
        if (IsIntrApproach)
        {
            Vector exetax(EMPTYEXTR); ObtainExtr(x, etax, &exetax);
            
            Vector tmp1(p, p, "complex"); tmp1.AlphaABaddBetaThis(GLOBAL::ZONE, x, GLOBAL::C, x, GLOBAL::N, GLOBAL::ZZERO);
            Vector tmp2(p, p, "complex"); tmp2.AlphaABaddBetaThis(GLOBAL::ZONE, x, GLOBAL::C, exetax, GLOBAL::N, GLOBAL::ZZERO);
            Vector S = tmp1 % tmp2; /* S = (x.GetTranspose() * x) % (x.GetTranspose() * exetax); */
            Vector Yp(exetax); Yp.AlphaABaddBetaThis(GLOBAL::ZNONE, x, GLOBAL::N, S, GLOBAL::N, GLOBAL::ZONE); /* Yp = etax - x * S; */
            
            Vector YYp(n, 2 * p, "complex");
            YYp.SubmatrixAssignment(0, n - 1, 0, p - 1, x);
            YYp.SubmatrixAssignment(0, n - 1, p, 2 * p - 1, Yp);
            YYp.QRDecom();
            Vector YYp_Q = YYp.Field("_Q"), YYp_R = YYp.Field("_R");

            Vector Z(2 * p, 2 * p, "complex"), Identity(p, p, "complex");
            Z.SetToZeros();
            Identity.SetToIdentity();
            Z.SubmatrixAssignment(0, p - 1, 0, p - 1, Identity + S + S.GetTranspose());
            Z.SubmatrixAssignment(0, p - 1, p, 2 * p - 1, Identity);
            Z.SubmatrixAssignment(p, 2 * p - 1, 0, p - 1, Identity);
            Z = YYp_R * Z * YYp_R.GetTranspose();
            Z.EigenDecomSym();
            Vector ZD = Z.Field("_EigVal").GetSubmatrix(p, 2 * p - 1, 0, 0).GetSqrt(), ZU = Z.Field("_EigVec").GetSubmatrix(0, 2 * p - 1, p, 2 * p - 1);
            Vector ZUD = ZD.GetDiagTimesM(ZU, GLOBAL::R);
            
            /*Find a representation Y such that Y minimize \|X - Y\|_F*/
            Yp.AlphaABaddBetaThis(GLOBAL::ZONE, YYp_Q, GLOBAL::N, ZUD, GLOBAL::N, GLOBAL::ZZERO); /* Yp = YYp_Q * ZUD; */
            Vector XtY(p, p, "complex"); XtY.AlphaABaddBetaThis(GLOBAL::ZONE, x, GLOBAL::C, Yp, GLOBAL::N, GLOBAL::ZZERO); /* XtY = x.GetTranspose() * Yp; */
            XtY.SVDDecom();
            tmp1.AlphaABaddBetaThis(GLOBAL::ZONE, XtY.Field("_Vt"), GLOBAL::C, XtY.Field("_U"), GLOBAL::C, GLOBAL::ZZERO);
            result->AlphaABaddBetaThis(GLOBAL::ZONE, Yp, GLOBAL::N, tmp1, GLOBAL::N, GLOBAL::ZZERO); /* result = Yp * XtY.Field("_Vt").GetTranspose() * XtY.Field("_U").GetTranspose(); */
            return *result;
        }
        
        return Manifold::Retraction(x, etax, result);
    };

    Vector &CSymFixedRankQ::DiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir) const
    {/* this is the differentiation of the retraction R_x(etax) = x + etax, not the one by projection */
        realdp nxix = std::sqrt(Metric(x, xix, xix));
        if (IsIntrApproach)
        {
            Vector exxix(EMPTYEXTR); ObtainExtr(x, xix, &exxix);
            ExtrProjection(y, exxix, &exxix);
            ObtainIntr(y, exxix, result);
        }
        else
        {
            ExtrProjection(y, xix, result);
        }

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

    Vector &CSymFixedRankQ::EucGradToGrad(const Variable &x, const Vector &egf, const Problem *prob, Vector *result)  const
    {
        if (metric == CSFRANKQEUC)
            return EucGradToGradEUC(x, egf, prob, result);

        if (metric == CSFRANKQHGZ)
            return EucGradToGradHGZ(x, egf, prob, result);
        
        printf("warning: CSymFixedRankQ::EucGradToGrad has not been done for this metric!\n");
        return Projection(x, egf, result);
    };

    Vector &CSymFixedRankQ::EucHvToHv(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const
    {
        if (metric == CSFRANKQEUC)
            return EucHvToHvEUC(x, etax, exix, prob, result);
        
        if (metric == CSFRANKQHGZ)
            return EucHvToHvHGZ(x, etax, exix, prob, result);
            
        printf("warning: CSymFixedRankQ::EucHvToHv has not been done for this metric!\n");
        return Projection(x, exix, result);
    };

    Vector &CSymFixedRankQ::EucGradToGradEUC(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const
    {
        if (prob->GetUseHess())
        {
            /*The copy on write is necessary. The reason is that the egf may be from a component in a product of elements.
            Therefore, if CopyOnWrite is not used, then the attached data in x and the product of elements share the same
            memory. This may cause an issue: if the product of elements are released before the attached data in x, then
            release the attached data in x would attempt to delete memory that has been released. This is an error!*/
            Vector Sharedegf(egf);
            Sharedegf.CopyOnWrite();
            x.AddToFields("EGrad", Sharedegf);
        }
        
        Vector XTX(p, p, "complex"); XTX.AlphaABaddBetaThis(GLOBAL::ZONE, x, GLOBAL::C, x, GLOBAL::N, GLOBAL::ZZERO);
        *result = (egf / XTX) / 2; /* result = (egf / (x.GetTranspose() * x)) / 2 */
        Vector tmp1(p, p, "complex"); tmp1.AlphaABaddBetaThis(GLOBAL::ZONE, x, GLOBAL::C, *result, GLOBAL::N, GLOBAL::ZZERO); /*tmp1 = x.GetTranspose() * ((egf / XTX) / 2) */
        tmp1 = XTX % tmp1;
        realdpcomplex coef = {-0.5, 0};
        result->AlphaABaddBetaThis(coef, x, GLOBAL::N, tmp1, GLOBAL::N, GLOBAL::ZONE);
        return *result;
    };

    /*for EUC metric*/
    Vector &CSymFixedRankQ::EucHvToHvEUC(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const
    {
        Vector egf = x.Field("EGrad");
        
        Vector XTX(p, p, "complex"); XTX.AlphaABaddBetaThis(GLOBAL::ZONE, x, GLOBAL::C, x, GLOBAL::N, GLOBAL::ZZERO);
        Vector EGYYinv = egf / XTX;
        Vector xTEGYYinv(p, p, "complex"); xTEGYYinv.AlphaABaddBetaThis(GLOBAL::ZONE, x, GLOBAL::C, EGYYinv, GLOBAL::N, GLOBAL::ZZERO);
        Vector YYinvYTEGYYinv = XTX % xTEGYYinv;
        
        Vector EtaTEGYYinv(p, p, "complex"); EtaTEGYYinv.AlphaABaddBetaThis(GLOBAL::ZONE, etax, GLOBAL::C, EGYYinv, GLOBAL::N, GLOBAL::ZZERO);
        Vector YTeta(p, p, "complex"); YTeta.AlphaABaddBetaThis(GLOBAL::ZONE, x, GLOBAL::C, etax, GLOBAL::N, GLOBAL::ZZERO);
        Vector YTxix(p, p, "complex"); YTxix.AlphaABaddBetaThis(GLOBAL::ZONE, x, GLOBAL::C, exix, GLOBAL::N, GLOBAL::ZZERO);
        Vector YYinvYTEH = XTX % YTxix;
        
        /*xix =  EucHess - 0.5 * Y * (Y^T * Y)^{-1} * Y^T * EucHess - 0.5 * Y * (Y^T * Y)^{-1} * Egrad^T * eta - Egrad * (Y^T * Y)^{-1} * Y * eta
                 + Y * (Y^T * Y)^{-1} * Y^T * Egrad * (Y^T * Y)^{-1} * Y^T * eta */
        realdpcomplex coef = {-0.5, 0};
        *result = exix;
        result->AlphaABaddBetaThis(coef, x, GLOBAL::N, YYinvYTEH, GLOBAL::N, GLOBAL::ZONE);
        result->AlphaABaddBetaThis(coef, x, GLOBAL::N, EtaTEGYYinv, GLOBAL::C, GLOBAL::ZONE);
        result->AlphaABaddBetaThis(GLOBAL::ZNONE, EGYYinv, GLOBAL::N, YTeta, GLOBAL::N, GLOBAL::ZONE);
        xTEGYYinv.AlphaABaddBetaThis(GLOBAL::ZONE, YYinvYTEGYYinv, GLOBAL::N, YTeta, GLOBAL::N, GLOBAL::ZZERO); /*xTEGYYinv is for temporary storage*/
        result->AlphaABaddBetaThis(GLOBAL::ZONE, x, GLOBAL::N, xTEGYYinv, GLOBAL::N, GLOBAL::ZONE); /* result  = exix - 0.5 * x * YYinvYTEH - 0.5 * x * EtaTEGYYinv.GetTranspose() - EGYYinv * YTeta + x * (YYinvYTEGYYinv * YTeta); */
        EGYYinv = 0.5 * ((*result) / XTX); /*EGYYinv is for temporary storage*/
        
        return ExtrProjection(x, EGYYinv, result);
    };

    /*for HGZ metric*/
    Vector &CSymFixedRankQ::EucGradToGradHGZ(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const
    {
        if (prob->GetUseHess())
        {
            /*The copy on write is necessary. The reason is that the egf may be from a component in a product of elements.
            Therefore, if CopyOnWrite is not used, then the attached data in x and the product of elements share the same
            memory. This may cause an issue: if the product of elements are released before the attached data in x, then
            release the attached data in x would attempt to delete memory that has been released. This is an error!*/
            Vector Sharedegf(egf);
            Sharedegf.CopyOnWrite();
            x.AddToFields("EGrad", Sharedegf);
        }
        Vector XTX(p, p, "complex"); XTX.AlphaABaddBetaThis(GLOBAL::ZONE, x, GLOBAL::C, x, GLOBAL::N, GLOBAL::ZZERO);
        *result = egf / XTX;
        return *result;
    };

    /*for HGZ metric*/
    Vector &CSymFixedRankQ::EucHvToHvHGZ(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const
    {
        Vector egf = x.Field("EGrad");
        
        Vector XTX(p, p, "complex"); XTX.AlphaABaddBetaThis(GLOBAL::ZONE, x, GLOBAL::C, x, GLOBAL::N, GLOBAL::ZZERO);
        Vector EGYYinv = egf / XTX;
        Vector YTEGYYinv(p, p, "complex"); YTEGYYinv.AlphaABaddBetaThis(GLOBAL::ZONE, x, GLOBAL::C, EGYYinv, GLOBAL::N, GLOBAL::ZZERO);
        Vector EtaTEGYYinv(p, p, "complex"); EtaTEGYYinv.AlphaABaddBetaThis(GLOBAL::ZONE, etax, GLOBAL::C, EGYYinv, GLOBAL::N, GLOBAL::ZZERO);
        Vector EtaTY(p, p, "complex"); EtaTY.AlphaABaddBetaThis(GLOBAL::ZONE, etax, GLOBAL::C, x, GLOBAL::N, GLOBAL::ZZERO);
        /*xix =  EucHess - 0.5 * Egrad * (Y^T * Y)^{-1} * Y^T * eta - 0.5 * Y * (Y^T * Y)^{-1} Egrad^T * eta + 0.5 * eta * Y^T * Egrad * (Y^T * Y)^{-1}
                 - 0.5 * Y * eta^T * Egrad * (Y^T * Y)^{-1} + 0.5 * eta *(Y^T * Y)^{-1} Egrad^T Y - 0.5 * Egrad * (Y^T * Y)^{-1} * eta^T * Y */
        *result = exix;
        realdpcomplex half = {0.5, 0}, nhalf = {-0.5, 0};
        result->AlphaABaddBetaThis(nhalf, EGYYinv, GLOBAL::N, EtaTY, GLOBAL::C, GLOBAL::ZONE);
        result->AlphaABaddBetaThis(nhalf, x, GLOBAL::N, EtaTEGYYinv, GLOBAL::C, GLOBAL::ZONE);
        result->AlphaABaddBetaThis(half, etax, GLOBAL::N, YTEGYYinv, GLOBAL::N, GLOBAL::ZONE);
        result->AlphaABaddBetaThis(nhalf, x, GLOBAL::N, EtaTEGYYinv, GLOBAL::N, GLOBAL::ZONE);
        result->AlphaABaddBetaThis(half, etax, GLOBAL::N, YTEGYYinv, GLOBAL::C, GLOBAL::ZONE);
        result->AlphaABaddBetaThis(nhalf, EGYYinv, GLOBAL::N, EtaTY, GLOBAL::N, GLOBAL::ZONE); /*result = exix - 0.5 * EGYYinv * EtaTY.GetTranspose() - 0.5 * x * EtaTEGYYinv.GetTranspose() + 0.5 * etax * YTEGYYinv - 0.5 * x * EtaTEGYYinv + 0.5 * etax * YTEGYYinv.GetTranspose() - 0.5 * EGYYinv * EtaTY;*/
        EGYYinv = *result / XTX; /*EGYYinv = result / (x.GetTranspose() * x); EGYYinv is used for temporaray storage*/
        return ExtrProjection(x, EGYYinv, result);
    };

    /*for EUC metric*/
    Vector &CSymFixedRankQ::ObtainIntrEUC(const Variable &x, const Vector &etax, Vector *result) const
    {
        if(!x.FieldsExist("_HHR"))
        {
            x.HHRDecom();
        }
        
        Vector tmp = etax.HHRMtp(x.Field("_HHR"), x.Field("_tau"), GLOBAL::C, GLOBAL::L);
        Vector HHR = x.Field("_HHR");
        const realdp *HHRptr = HHR.ObtainReadData();
        realdp *tmpptr = tmp.ObtainWritePartialData();
        for (integer i = 0; i < p; i++)
        {
            if(((realdpcomplex *) HHRptr)[i + n * i].r < 0)
                scal_(&p, &GLOBAL::ZNONE, ((realdpcomplex *) tmpptr) + i, &n);
        }
        Vector L(p, p, "complex");
        realdpcomplex *Lptr = (realdpcomplex *) L.ObtainWriteEntireData();
        for(integer i = 0; i < p; i++)
        {
            for(integer j = 0; j < i; j++)
            {
                Lptr[j + i * p].r = 0;
                Lptr[j + i * p].i = 0;
            }
            for(integer j = i; j < p; j++)
            {
                Lptr[j + i * p].r = ((realdpcomplex *) HHRptr)[i + n * j].r;
                Lptr[j + i * p].i = - (((realdpcomplex *) HHRptr)[i + n * j].i);
                if(((realdpcomplex *)HHRptr)[i + n * i].r < 0)
                {
                    Lptr[j + i * p].r *= -1;
                    Lptr[j + i * p].i *= -1;
                }
            }
        }
        tmp = tmp * L;
        tmpptr = tmp.ObtainWritePartialData(); /*tmp has been updated. we need to reload its pointer*/
        realdp *resultptr = result->ObtainWriteEntireData();
        integer idx = 0;
        realdp r2 = static_cast<realdp> (sqrt(2.0));
        
        for (integer i = 0; i < p; i++)
        {
            resultptr[idx] = 2.0 * ((realdpcomplex *) tmpptr)[i + i * n].r;
            idx++;
        }
        for (integer i = 0; i < p; i++)
        {
            for (integer j = i + 1; j < p; j++)
            {
                resultptr[idx] = 2.0 * r2 * ((realdpcomplex *) tmpptr)[j + i * n].r;
                idx++;
                resultptr[idx] = 2.0 * r2 * ((realdpcomplex *) tmpptr)[j + i * n].i;
                idx++;
            }
        }

        for (integer i = 0; i < p; i++)
        {
            for (integer j = p; j < n; j++)
            {
                resultptr[idx] = r2 * ((realdpcomplex *) tmpptr)[j + i * n].r;
                idx++;
                resultptr[idx] = r2 * ((realdpcomplex *) tmpptr)[j + i * n].i;
                idx++;
            }
        }
        return *result;
    };

    /*for EUC metric*/
    Vector &CSymFixedRankQ::ObtainExtrEUC(const Variable &x, const Vector &intretax, Vector *result) const
    {
        if(!x.FieldsExist("_HHR"))
        {
            x.HHRDecom();
        }
        Vector HHR = x.Field("_HHR");
        const realdp *HHRptr = HHR.ObtainReadData();
        realdp *resultptr = result->ObtainWriteEntireData();
        const realdp *intretaxptr = intretax.ObtainReadData();
        realdp r2 = static_cast<realdp> (sqrt(2.0));
        integer idx = 0;
        for (integer i = 0; i < p; i++)
        {
            ((realdpcomplex *)resultptr)[i + i * n].r = intretaxptr[idx] / 2.0;
            ((realdpcomplex *)resultptr)[i + i * n].i = 0;
            idx++;
        }
        
        for (integer i = 0; i < p; i++)
        {
            for (integer j = i + 1; j < p; j++)
            {
                ((realdpcomplex *)resultptr)[j + i * n].r = intretaxptr[idx] / r2 / 2.0;
                ((realdpcomplex *)resultptr)[i + j * n].r = intretaxptr[idx] / r2 / 2.0;
                idx++;
                ((realdpcomplex *)resultptr)[j + i * n].i = intretaxptr[idx] / r2 / 2.0;
                ((realdpcomplex *)resultptr)[i + j * n].i = -intretaxptr[idx] / r2 / 2.0;
                idx++;
            }
        }
        for (integer i = 0; i < p; i++)
        {
            for (integer j = p; j < n; j++)
            {
                ((realdpcomplex *)resultptr)[j + i * n].r = intretaxptr[idx] / r2;
                idx++;
                ((realdpcomplex *)resultptr)[j + i * n].i = intretaxptr[idx] / r2;
                idx++;
            }
        }
        Vector L(p, p, "complex");
        realdpcomplex *Lptr = (realdpcomplex *) L.ObtainWriteEntireData();
        for(integer i = 0; i < p; i++)
        {
            for(integer j = 0; j < i; j++)
            {
                Lptr[j + i * p].r = 0;
                Lptr[j + i * p].i = 0;
            }
            for(integer j = i; j < p; j++)
            {
                Lptr[j + i * p].r = ((realdpcomplex *) HHRptr)[i + n * j].r;
                Lptr[j + i * p].i = - (((realdpcomplex *) HHRptr)[i + n * j].i);
                if(((realdpcomplex *)HHRptr)[i + n * i].r < 0)
                {
                    Lptr[j + i * p].r *= -1;
                    Lptr[j + i * p].i *= -1;
                }
            }
        }
        *result = *result / L;
        resultptr = result->ObtainWritePartialData();
        
        for (integer i = 0; i < p; i++)
        {
            if(((realdpcomplex *)HHRptr)[i + n * i].r < 0)
                scal_(&p, &GLOBAL::ZNONE, ((realdpcomplex *)resultptr) + i, &n);
        }
        *result = result->HHRMtp(x.Field("_HHR"), x.Field("_tau"), GLOBAL::N, GLOBAL::L);
        return *result;
    };


    /*for HGZ metric*/
    Vector &CSymFixedRankQ::ObtainIntrHGZ(const Variable &x, const Vector &etax, Vector *result) const
    {
        if(!x.FieldsExist("_HHR"))
        {
            x.HHRDecom();
        }
        
        Vector tmp = etax.HHRMtp(x.Field("_HHR"), x.Field("_tau"), GLOBAL::C, GLOBAL::L);
        Vector HHR = x.Field("_HHR");
        const realdp *HHRptr = HHR.ObtainReadData();
        realdp *tmpptr = tmp.ObtainWritePartialData();
        for (integer i = 0; i < p; i++)
        {
            if(((realdpcomplex *) HHRptr)[i + n * i].r < 0)
                scal_(&p, &GLOBAL::ZNONE, ((realdpcomplex *) tmpptr) + i, &n);
        }
        Vector L(p, p, "complex");
        realdpcomplex *Lptr = (realdpcomplex *) L.ObtainWriteEntireData();
        for(integer i = 0; i < p; i++)
        {
            for(integer j = 0; j < i; j++)
            {
                Lptr[j + i * p].r = 0;
                Lptr[j + i * p].i = 0;
            }
            for(integer j = i; j < p; j++)
            {
                Lptr[j + i * p].r = ((realdpcomplex *) HHRptr)[i + n * j].r;
                Lptr[j + i * p].i = - (((realdpcomplex *) HHRptr)[i + n * j].i);
                if(((realdpcomplex *)HHRptr)[i + n * i].r < 0)
                {
                    Lptr[j + i * p].r *= -1;
                    Lptr[j + i * p].i *= -1;
                }
            }
        }
        
        tmp = tmp * L;
        tmpptr = tmp.ObtainWritePartialData(); /*tmp has been updated. we need to reload its pointer*/
        realdp *resultptr = result->ObtainWriteEntireData();
        integer idx = 0;
        realdp r2 = static_cast<realdp> (sqrt(2.0));
        
        for (integer i = 0; i < p; i++)
        {
            resultptr[idx] = ((realdpcomplex *) tmpptr)[i + i * n].r;
            idx++;
        }
        for (integer i = 0; i < p; i++)
        {
            for (integer j = i + 1; j < p; j++)
            {
                resultptr[idx] = r2 * ((realdpcomplex *) tmpptr)[j + i * n].r;
                idx++;
                resultptr[idx] = r2 * ((realdpcomplex *) tmpptr)[j + i * n].i;
                idx++;
            }
        }

        for (integer i = 0; i < p; i++)
        {
            for (integer j = p; j < n; j++)
            {
                resultptr[idx] = ((realdpcomplex *) tmpptr)[j + i * n].r;
                idx++;
                resultptr[idx] = ((realdpcomplex *) tmpptr)[j + i * n].i;
                idx++;
            }
        }
        return *result;
    };

    /*for HGZ metric*/
    Vector &CSymFixedRankQ::ObtainExtrHGZ(const Variable &x, const Vector &intretax, Vector *result) const
    {
        if(!x.FieldsExist("_HHR"))
        {
            x.HHRDecom();
        }
        Vector HHR = x.Field("_HHR");
        const realdp *HHRptr = HHR.ObtainReadData();
        realdp *resultptr = result->ObtainWriteEntireData();
        const realdp *intretaxptr = intretax.ObtainReadData();
        realdp r2 = static_cast<realdp> (sqrt(2.0));
        integer idx = 0;
        for (integer i = 0; i < p; i++)
        {
            ((realdpcomplex *)resultptr)[i + i * n].r = intretaxptr[idx];
            ((realdpcomplex *)resultptr)[i + i * n].i = 0;
            idx++;
        }

        for (integer i = 0; i < p; i++)
        {
            for (integer j = i + 1; j < p; j++)
            {
                ((realdpcomplex *)resultptr)[j + i * n].r = intretaxptr[idx] / r2;
                ((realdpcomplex *)resultptr)[i + j * n].r = intretaxptr[idx] / r2;
                idx++;
                ((realdpcomplex *)resultptr)[j + i * n].i = intretaxptr[idx] / r2;
                ((realdpcomplex *)resultptr)[i + j * n].i = -intretaxptr[idx] / r2;
                idx++;
            }
        }
        for (integer i = 0; i < p; i++)
        {
            for (integer j = p; j < n; j++)
            {
                ((realdpcomplex *)resultptr)[j + i * n].r = intretaxptr[idx];
                idx++;
                ((realdpcomplex *)resultptr)[j + i * n].i = intretaxptr[idx];
                idx++;
            }
        }
        Vector L(p, p, "complex");
        realdpcomplex *Lptr = (realdpcomplex *) L.ObtainWriteEntireData();
        for(integer i = 0; i < p; i++)
        {
            for(integer j = 0; j < i; j++)
            {
                Lptr[j + i * p].r = 0;
                Lptr[j + i * p].i = 0;
            }
            for(integer j = i; j < p; j++)
            {
                Lptr[j + i * p].r = ((realdpcomplex *) HHRptr)[i + n * j].r;
                Lptr[j + i * p].i = - (((realdpcomplex *) HHRptr)[i + n * j].i);
                if(((realdpcomplex *)HHRptr)[i + n * i].r < 0)
                {
                    Lptr[j + i * p].r *= -1;
                    Lptr[j + i * p].i *= -1;
                }
            }
        }
        *result = *result / L;
        resultptr = result->ObtainWritePartialData();
        
        for (integer i = 0; i < p; i++)
        {
            if(((realdpcomplex *)HHRptr)[i + n * i].r < 0)
                scal_(&p, &GLOBAL::ZNONE, ((realdpcomplex *)resultptr) + i, &n);
        }
        
        *result = result->HHRMtp(x.Field("_HHR"), x.Field("_tau"), GLOBAL::N, GLOBAL::L);
        return *result;
    };

    /*The vector transport by projection*/
    Vector &CSymFixedRankQ::VectorTransportProj(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const
    {
        if(IsIntrApproach)
        {
            Vector exxix(EMPTYEXTR); ObtainExtr(x, xix, &exxix);
            return ObtainIntr(y, exxix, result);
        }
        
        return ExtrProjection(y, xix, result);
    };
}; /*end of ROPTLIB namespace*/
