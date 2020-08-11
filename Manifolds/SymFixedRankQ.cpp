
#include "Manifolds/SymFixedRankQ.h"

/*Define the namespace*/
namespace ROPTLIB{

	SymFixedRankQ::SymFixedRankQ(integer r, integer c)
	{
		n = r;
		p = c;
		IsIntrApproach = true;
		HasHHR = false;
		name.assign("Symmetric fixed rank manifold by quotient representation");
		IntrinsicDim = r * c - c * (c - 1) / 2;
		ExtrinsicDim = r * c;
        metric = SFRANKQEUC;
        VecTran = SFRANKQPARALLELIZATION;
        EMPTYEXTR = Vector(r, c);
		EMPTYINTR = Vector(IntrinsicDim);
	};

	void SymFixedRankQ::SetParams(PARAMSMAP params)
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
				case 5:
					ChooseParamsSet5();
					break;
				case 6:
					ChooseParamsSet6();
					break;
				default:
					printf("warning: invalid index for parameter set! use the default one.");
					ChooseParamsSet1();
					break;
				}
			}
		}
	};

	void SymFixedRankQ::ChooseParamsSet1(void)
	{
        metric = SFRANKQEUC;
        VecTran = SFRANKQPARALLELIZATION;
	};

	void SymFixedRankQ::ChooseParamsSet2(void)
	{
        metric = SFRANKQHGZ;
        VecTran = SFRANKQPARALLELIZATION;
	};

	void SymFixedRankQ::ChooseParamsSet3(void)
	{
		metric = SFRANKQQEUC;
		VecTran = SFRANKQPARALLELIZATION;
	};

	void SymFixedRankQ::ChooseParamsSet4(void)
	{
        metric = SFRANKQEUC;
        VecTran = SFRANKQPROJECTION;
	};

	void SymFixedRankQ::ChooseParamsSet5(void)
	{
        metric = SFRANKQHGZ;
        VecTran = SFRANKQPROJECTION;
	};

	void SymFixedRankQ::ChooseParamsSet6(void)
	{
		metric = SFRANKQQEUC;
		VecTran = SFRANKQPROJECTION;
	};

	SymFixedRankQ::~SymFixedRankQ(void)
	{
	};

    Variable SymFixedRankQ::RandominManifold(void) const
    {
        Variable result(EMPTYEXTR); result.RandGaussian();
        return result;
    };

    void SymFixedRankQ::CheckParams(void) const
    {
        std::string SFRANKQMetricnames[SFRANKQMETRICLENGTH] = { "SFRANKQEUC", "SFRANKQQEUC", "SFRANKQHGZ" };
        std::string SFRANKQVTnames[SFRANKQVECTORTRANSPORTLENGTH] = { "SFRANKQPARALLELIZATION", "SFRANKQPROJECTION" };
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

	realdp SymFixedRankQ::Metric(const Variable &x, const Vector &etax, const Vector &xix) const
	{
		if (IsIntrApproach)
			return Manifold::Metric(x, etax, xix);
        
        if(metric == SFRANKQEUC)
        {/*2 \trace(x^T etax x^T xix + x^T x etax^T xix)*/
            Vector tmp1(p, p); tmp1.AlphaABaddBetaThis(1, etax, GLOBAL::T, x, GLOBAL::N, 0);
            Vector tmp2(p, p); tmp2.AlphaABaddBetaThis(1, x, GLOBAL::T, xix, GLOBAL::N, 0);
            realdp result = tmp1.DotProduct(tmp2);
            tmp1.AlphaABaddBetaThis(1, x, GLOBAL::T, x, GLOBAL::N, 0);
            tmp2.AlphaABaddBetaThis(1, etax, GLOBAL::T, xix, GLOBAL::N, 0);
            result += tmp1.DotProduct(tmp2);
            return 2.0 * result;
        }
        
        if(metric == SFRANKQQEUC)
        {/*\trace(etax^T xix)*/
            return Manifold::Metric(x, etax, xix);
        }
        
        if(metric == SFRANKQHGZ)
        { /*\trace(x^T x etax^T xix)*/
            Vector tmp1(p, p); tmp1.AlphaABaddBetaThis(1, x, GLOBAL::T, x, GLOBAL::N, 0);
            Vector tmp2(p, p); tmp2.AlphaABaddBetaThis(1, etax, GLOBAL::T, xix, GLOBAL::N, 0);
            return tmp1.DotProduct(tmp2);
        }
        
		printf("Warning: this metric has not been done!\n");
		return 0;
	};

    Vector &SymFixedRankQ::Projection(const Variable &x, const Vector &etax, Vector *result) const
    {
        if (IsIntrApproach)
            return IntrProjection(x, etax, result);
        
        return ExtrProjection(x, etax, result);
    };

    Vector &SymFixedRankQ::ExtrProjection(const Variable &x, const Vector &etax, Vector *result) const
    {
        if (metric == SFRANKQEUC || metric == SFRANKQHGZ)
            return ExtrProjectionEUCorHGZ(x, etax, result);
        
        if (metric == SFRANKQQEUC)
            return ExtrProjectionQEUC(x, etax, result);
        
        printf("Warning: SymFixedRankQ::ExtrProjection has not been done for this Riemannian metric!\n");
        return Manifold::ExtrProjection(x, etax, result);
    };

    Vector &SymFixedRankQ::ExtrProjectionQEUC(const Variable &x, const Vector &etax, Vector *result) const
    {
        Vector XTX(p, p); XTX.AlphaABaddBetaThis(1, x, GLOBAL::T, x, GLOBAL::N, 0);
        Vector XTV(p, p); XTV.AlphaABaddBetaThis(1, x, GLOBAL::T, etax, GLOBAL::N, 0);
        XTV = XTV - XTV.GetTranspose();
        
        Vector tmp(XTV.SYL(XTX, XTX));
        *result = etax;
        result->AlphaABaddBetaThis(-1, x, GLOBAL::N, tmp, GLOBAL::N, 1); /* result = etax - x * XTV.SYL(XTX, XTX)*/
        return *result;
    };

    Vector &SymFixedRankQ::ExtrProjectionEUCorHGZ(const Variable &x, const Vector &etax, Vector *result) const
    {
        Vector XTX(p, p); XTX.AlphaABaddBetaThis(1, x, GLOBAL::T, x, GLOBAL::N, 0);
        Vector XTV(p, p); XTV.AlphaABaddBetaThis(1, x, GLOBAL::T, etax, GLOBAL::N, 0);
        XTV = XTX % XTV;
        
        *result = etax;
        XTX = (XTV - XTV.GetTranspose()) / 2.0;
        result->AlphaABaddBetaThis(-1, x, GLOBAL::N, XTX, GLOBAL::N, 1); /* result = etax - x * ((XTV - XTV.GetTranspose()) / 2.0);*/
        return *result;
    };

    Vector &SymFixedRankQ::ObtainIntr(const Variable &x, const Vector &etax, Vector *result) const
    {
        if (metric == SFRANKQEUC)
            return ObtainIntrEUC(x, etax, result);
        
        if (metric == SFRANKQHGZ)
            return ObtainIntrHGZ(x, etax, result);
            
        if (metric == SFRANKQQEUC)
            return ObtainIntrQEUC(x, etax, result);
        
        printf("warning: SymFixedRankQ::ObtainIntr has not been done for this metric!\n");
        return Manifold::ObtainIntr(x, etax, result);
    };

    Vector &SymFixedRankQ::ObtainExtr(const Variable &x, const Vector &intretax, Vector *result) const
    {
        if (metric == SFRANKQEUC)
            return ObtainExtrEUC(x, intretax, result);
        
        if (metric == SFRANKQHGZ)
            return ObtainExtrHGZ(x, intretax, result);
        
        if (metric == SFRANKQQEUC)
            return ObtainExtrQEUC(x, intretax, result);
        printf("warning: SymFixedRankQ::ObtainExtr has not been done for this metric!\n");
        return Manifold::ObtainExtr(x, intretax, result);
    };

    Vector &SymFixedRankQ::coTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
    {
        printf("warning:SymFixedRankQ::coTangentVector has not been done!\n");
        return Manifold::coTangentVector(x, etax, y, xiy, result);
    };

    Vector &SymFixedRankQ::VectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const
    {
        if (VecTran == SFRANKQPARALLELIZATION && !HasHHR)
        {
            return Manifold::VectorTransport(x, etax, y, xix, result);
        }

        if (VecTran == SFRANKQPROJECTION && !HasHHR)
        {
            return VectorTransportProj(x, etax, y, xix, result);
        }

        if (HasHHR)
            return LCVectorTransport(x, etax, y, xix, result);

        printf("Error: VectorTransport has not been done!\n");
        return Manifold::VectorTransport(x, etax, y, xix, result);
    };

    Vector &SymFixedRankQ::InverseVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result)  const
    {
        if (VecTran == SFRANKQPARALLELIZATION && !HasHHR)
        {
            return Manifold::InverseVectorTransport(x, etax, y, xiy, result);
        }

        if (VecTran == SFRANKQPROJECTION && !HasHHR)
        {
            printf("Warning: Inverse Vector transport by projection has not been done!\n");
            return Manifold::InverseVectorTransport(x, etax, y, xiy, result);
        }

        if (HasHHR)
            return LCInverseVectorTransport(x, etax, y, xiy, result);

        printf("Error: InverseVectorTransport has not been done!\n");
        return Manifold::InverseVectorTransport(x, etax, y, xiy, result);
    };

    LinearOPE &SymFixedRankQ::TranHInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, LinearOPE *result) const
    {
        if (VecTran == SFRANKQPARALLELIZATION && !HasHHR)
        {
            return Manifold::TranHInvTran(x, etax, y, Hx, result);
        }

        if (VecTran == SFRANKQPROJECTION && !HasHHR)
        {
            printf("Warning: Transport a linear operator using vector transport by projection has not been done!\n");
            return Manifold::TranHInvTran(x, etax, y, Hx, result);
        }

        if (HasHHR)
            return LCVectorTransport(x, etax, y, Hx, result);

        printf("Error: TranHInvTran has not been done!\n");
        return Manifold::TranHInvTran(x, etax, y, Hx, result);
    };

    Variable &SymFixedRankQ::Retraction(const Variable &x, const Vector &etax, Variable *result) const
    {
        if (IsIntrApproach)
        {
            Vector exetax(EMPTYEXTR);
            ObtainExtr(x, etax, &exetax);

            Vector tmp1(p, p); tmp1.AlphaABaddBetaThis(1, x, GLOBAL::T, x, GLOBAL::N, 0);
            Vector tmp2(p, p); tmp2.AlphaABaddBetaThis(1, x, GLOBAL::T, exetax, GLOBAL::N, 0);
            Vector S = tmp1 % tmp2; /* S = (x.GetTranspose() * x) % (x.GetTranspose() * exetax); */
            Vector Yp(exetax); Yp.AlphaABaddBetaThis(-1, x, GLOBAL::N, S, GLOBAL::N, 1); /* Yp = etax - x * S; */
            
            Vector YYp(n, 2 * p);
            YYp.SubmatrixAssignment(0, n - 1, 0, p - 1, x);
            YYp.SubmatrixAssignment(0, n - 1, p, 2 * p - 1, Yp);
            YYp.QRDecom();
            Vector YYp_Q = YYp.Field("_Q"), YYp_R = YYp.Field("_R");

            Vector Z(2 * p, 2 * p), Identity(p, p);
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
            Yp.AlphaABaddBetaThis(1, YYp_Q, GLOBAL::N, ZUD, GLOBAL::N, 0); /* Yp = YYp_Q * ZUD; */
            Vector XtY(p, p); XtY.AlphaABaddBetaThis(1, x, GLOBAL::T, Yp, GLOBAL::N, 0); /* XtY = x.GetTranspose() * Yp; */
            XtY.SVDDecom();
            tmp1.AlphaABaddBetaThis(1, XtY.Field("_Vt"), GLOBAL::T, XtY.Field("_U"), GLOBAL::T, 0);
            result->AlphaABaddBetaThis(1, Yp, GLOBAL::N, tmp1, GLOBAL::N, 0); /* result = Yp * XtY.Field("_Vt").GetTranspose() * XtY.Field("_U").GetTranspose(); */
            return *result;
        }
        
        return Manifold::Retraction(x, etax, result);
    };

    Vector &SymFixedRankQ::DiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir) const
    { /* this is the differentiation of the retraction R_x(etax) = x + etax, not the one by projection */
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

    Vector &SymFixedRankQ::EucGradToGrad(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const
	{
		if (metric == SFRANKQEUC)
			return EucGradToGradEUC(x, egf, prob, result);

        if (metric == SFRANKQHGZ)
            return EucGradToGradHGZ(x, egf, prob, result);
        
        if (metric == SFRANKQQEUC)
            return EucGradToGradQEUC(x, egf, prob, result);
			
        printf("warning: SymFixedRankQ::EucGradToGrad has not been done for this metric!\n");
        return Projection(x, egf, result);
	};

	Vector &SymFixedRankQ::EucHvToHv(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const
	{
		if (metric == SFRANKQEUC)
			return EucHvToHvEUC(x, etax, exix, prob, result);
		
        if (metric == SFRANKQHGZ)
            return EucHvToHvHGZ(x, etax, exix, prob, result);
			
        if (metric == SFRANKQQEUC)
            return EucHvToHvQEUC(x, etax, exix, prob, result);
        
        printf("warning: SymFixedRankQ::EucHvToHv has not been done for this metric!\n");
        return Projection(x, exix, result);
	};

	Vector &SymFixedRankQ::EucGradToGradEUC(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const
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
        
        Vector XTX(p, p); XTX.AlphaABaddBetaThis(1, x, GLOBAL::T, x, GLOBAL::N, 0);
        *result = (egf / XTX) / 2; /* result = (egf / (x.GetTranspose() * x)) / 2 */
        Vector tmp1(p, p); tmp1.AlphaABaddBetaThis(1, x, GLOBAL::T, *result, GLOBAL::N, 0); /*tmp1 = x.GetTranspose() * tmp1 */
        tmp1 = XTX % tmp1;
        result->AlphaABaddBetaThis(-0.5, x, GLOBAL::N, tmp1, GLOBAL::N, 1);
        return *result;
	};

	Vector &SymFixedRankQ::EucHvToHvEUC(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const
	{
        Vector egf = x.Field("EGrad");
        
        Vector XTX(p, p); XTX.AlphaABaddBetaThis(1, x, GLOBAL::T, x, GLOBAL::N, 0);
        Vector EGYYinv = egf / XTX;
        Vector xTEGYYinv(p, p); xTEGYYinv.AlphaABaddBetaThis(1, x, GLOBAL::T, EGYYinv, GLOBAL::N, 0);
        Vector YYinvYTEGYYinv = XTX % xTEGYYinv;
        
        Vector EtaTEGYYinv(p, p); EtaTEGYYinv.AlphaABaddBetaThis(1, etax, GLOBAL::T, EGYYinv, GLOBAL::N, 0);
        Vector YTeta(p, p); YTeta.AlphaABaddBetaThis(1, x, GLOBAL::T, etax, GLOBAL::N, 0);
        Vector YTxix(p, p); YTxix.AlphaABaddBetaThis(1, x, GLOBAL::T, exix, GLOBAL::N, 0);
        Vector YYinvYTEH = XTX % YTxix;
        
        /*xix =  EucHess - 0.5 * Y * (Y^T * Y)^{-1} * Y^T * EucHess - 0.5 * Y * (Y^T * Y)^{-1} * Egrad^T * eta - Egrad * (Y^T * Y)^{-1} * Y * eta
                 + Y * (Y^T * Y)^{-1} * Y^T * Egrad * (Y^T * Y)^{-1} * Y^T * eta */
        *result = exix;
        result->AlphaABaddBetaThis(-0.5, x, GLOBAL::N, YYinvYTEH, GLOBAL::N, 1);
        result->AlphaABaddBetaThis(-0.5, x, GLOBAL::N, EtaTEGYYinv, GLOBAL::T, 1);
        result->AlphaABaddBetaThis(-1, EGYYinv, GLOBAL::N, YTeta, GLOBAL::N, 1);
        xTEGYYinv.AlphaABaddBetaThis(1, YYinvYTEGYYinv, GLOBAL::N, YTeta, GLOBAL::N, 0); /*xTEGYYinv is for temporary storage*/
        result->AlphaABaddBetaThis(1, x, GLOBAL::N, xTEGYYinv, GLOBAL::N, 1); /* result  = exix - 0.5 * x * YYinvYTEH - 0.5 * x * EtaTEGYYinv.GetTranspose() - EGYYinv * YTeta + x * (YYinvYTEGYYinv * YTeta); */
        EGYYinv = 0.5 * (*result) / XTX; /*EGYYinv is for temporary storage*/
        
        return ExtrProjection(x, EGYYinv, result);
    };

	Vector &SymFixedRankQ::ObtainIntrEUC(const Variable &x, const Vector &etax, Vector *result) const
	{
        if(!x.FieldsExist("_HHR"))
        {
            x.HHRDecom();
        }
        
        Vector tmp = etax.HHRMtp(x.Field("_HHR"), x.Field("_tau"), GLOBAL::T, GLOBAL::L);
        Vector HHR = x.Field("_HHR");
        const realdp *HHRptr = HHR.ObtainReadData();
        realdp *tmpptr = tmp.ObtainWritePartialData();
        for (integer i = 0; i < p; i++)
        {
            if(HHRptr[i + n * i] < 0)
                scal_(&p, &GLOBAL::DNONE, tmpptr + i, &n);
        }
        Vector L(p, p);
        realdp *Lptr = L.ObtainWriteEntireData();
        for(integer i = 0; i < p; i++)
        {
            for(integer j = 0; j < i; j++)
            {
                Lptr[j + i * p] = 0;
            }
            for(integer j = i; j < p; j++)
            {
                Lptr[j + i * p] = HHRptr[i + n * j];
                if(HHRptr[i + n * i] < 0)
                    Lptr[j + i * p] *= -1;
            }
        }
        
        tmp = tmp * L;
        tmpptr = tmp.ObtainWritePartialData();
        realdp *resultptr = result->ObtainWriteEntireData();
        
        integer idx = 0;
        realdp r2 = static_cast<realdp> (sqrt(2.0));
        
        for (integer i = 0; i < p; i++)
        {
            resultptr[idx] = 2.0 * tmpptr[i + i * n];
            idx++;
        }
        for (integer i = 0; i < p; i++)
        {
            for (integer j = i + 1; j < p; j++)
            {
                resultptr[idx] = 2.0 * r2 * tmpptr[j + i * n];
                idx++;
            }
        }

        for (integer i = 0; i < p; i++)
        {
            for (integer j = p; j < n; j++)
            {
                resultptr[idx] = r2 * tmpptr[j + i * n];
                idx++;
            }
        }
        return *result;
	};

	Vector &SymFixedRankQ::ObtainExtrEUC(const Variable &x, const Vector &intretax, Vector *result) const
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
            resultptr[i + i * n] = intretaxptr[idx] / 2.0;
            idx++;
        }

        for (integer i = 0; i < p; i++)
        {
            for (integer j = i + 1; j < p; j++)
            {
                resultptr[j + i * n] = intretaxptr[idx] / r2 / 2.0;
                resultptr[i + j * n] = intretaxptr[idx] / r2 / 2.0;
                idx++;
            }
        }
        for (integer i = 0; i < p; i++)
        {
            for (integer j = p; j < n; j++)
            {
                resultptr[j + i * n] = intretaxptr[idx] / r2;
                idx++;
            }
        }
        
        for (integer i = 0; i < p; i++)
        {
            if(HHRptr[i + n * i] < 0)
                scal_(&p, &GLOBAL::DNONE, resultptr + i, &n);
        }
        Vector L(p, p);
        realdp *Lptr = L.ObtainWriteEntireData();
        for(integer i = 0; i < p; i++)
        {
            for(integer j = 0; j < i; j++)
            {
                Lptr[j + i * p] = 0;
            }
            for(integer j = i; j < p; j++)
            {
                Lptr[j + i * p] = HHRptr[i + n * j];
                if(HHRptr[i + n * i] < 0)
                    Lptr[j + i * p] *= -1;
            }
        }
        
        *result = *result / L;
        *result = result->HHRMtp(x.Field("_HHR"), x.Field("_tau"), GLOBAL::N, GLOBAL::L);
        return *result;
	};

	Vector &SymFixedRankQ::EucGradToGradQEUC(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const
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
		return ExtrProjection(x, egf, result);
	};

	Vector &SymFixedRankQ::EucHvToHvQEUC(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const
	{
		return ExtrProjection(x, exix, result);
	};

	Vector &SymFixedRankQ::ObtainIntrQEUC(const Variable &x, const Vector &etax, Vector *result) const
	{
        if(!x.FieldsExist("_HHR"))
        {
            x.HHRDecom();
        }
        
        Vector tmp = etax.HHRMtp(x.Field("_HHR"), x.Field("_tau"), GLOBAL::T, GLOBAL::L);
        Vector HHR = x.Field("_HHR");
        const realdp *HHRptr = HHR.ObtainReadData();
        realdp *tmpptr = tmp.ObtainWritePartialData();
        for (integer i = 0; i < p; i++)
        {
            if(HHRptr[i + n * i] < 0)
                scal_(&p, &GLOBAL::DNONE, tmpptr + i, &n);
        }
        
        Vector L(p, p);
        realdp *Lptr = L.ObtainWriteEntireData();
        for(integer i = 0; i < p; i++)
        {
            for(integer j = 0; j < i; j++)
            {
                Lptr[j + i * p] = 0;
            }
            for(integer j = i; j < p; j++)
            {
                Lptr[j + i * p] = HHRptr[i + n * j];
                if(HHRptr[i + n * i] < 0)
                    Lptr[j + i * p] *= -1;
            }
        }
        Vector Ltmp = L * tmp.GetSubmatrix(0, p - 1, 0, p - 1);
        const realdp *Ltmpptr = Ltmp.ObtainReadData();
        realdp *resultptr = result->ObtainWriteEntireData();
        integer idx = 0;
        for (integer i = 0; i < p; i++)
        {
            for (integer j = i; j < p; j++)
            {
                resultptr[idx] = Ltmpptr[j + i * p];
                idx++;
            }
        }
        for (integer i = 0; i < p; i++)
        {
            for (integer j = p; j < n; j++)
            {
                resultptr[idx] = tmpptr[j + i * n];
                idx++;
            }
        }

        /* Above: result is then R^{-T} S + K. Next, we first find the coefficients of the etax under a non-orthonormal basis.
        Then orthonormalize the basis and update the coefficient. To the end, we find the non-orthogonal basis and store it in B.
        Then use qr decomposition to find the matrix R. Then use the matrix R to find the coefficient vector under the orthonormalized
        basis. total complexity is O(p^6) + O(n p)*/
        
        Vector B(p * p, p * (p + 1) / 2);
        B.SetToZeros();
        realdp *Bptr = B.ObtainWritePartialData();
        integer info;
        idx = 0;
        for (integer i = 0; i < p; i++)
        {
            for (integer j = i; j < p; j++)
            {
                Bptr[j + i * p + idx * p * p] = 1;
                Bptr[i + j * p + idx * p * p] = 1;
                trtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &p, &p, Lptr, &p, Bptr + idx * p * p, &p, &info);
                idx++;
            }
        }
        
        B.HHRDecom();
        Vector BHHR = B.Field("_HHR");
        const realdp *BHHRptr = BHHR.ObtainReadData();
        
        integer dim = p * (p + 1) / 2;
        Vector BL(dim, dim);
        realdp sign;
        realdp *BLptr = BL.ObtainWriteEntireData();
        for (integer i = 0; i < dim; i++)
        {
            for (integer j = 0; j < i; j++)
            {
                BLptr[j + i * dim] = 0;
            }
            sign = static_cast<realdp> ((BHHRptr[i + p * p * i] >= 0) ? 1 : -1);
            for (integer j = i; j < dim; j++)
            {
                BLptr[j + i * dim] = BHHRptr[i + p * p * j] * sign;
            }
        }

        result->SubmatrixAssignment(0, dim - 1, 0, 0, BL.GetTranspose() * result->GetSubmatrix(0, dim - 1, 0, 0));
        
        x.AddToFields("BL", BL);
        return *result;
	};

	Vector &SymFixedRankQ::ObtainExtrQEUC(const Variable &x, const Vector &intretax, Vector *result) const
	{
        if(!x.FieldsExist("_HHR"))
        {
            x.HHRDecom();
        }
        
        integer dim = p * (p + 1) / 2;
        Vector BL = x.Field("BL");
        
        Vector tmp = BL.GetTranspose() % intretax.GetSubmatrix(0, dim - 1, 0, 0);
        const realdp *tmpptr = tmp.ObtainReadData();
        const realdp *intretaxptr = intretax.ObtainReadData();
        
        realdp *resultptr = result->ObtainWriteEntireData();
        integer idx = 0;
        for (integer i = 0; i < p; i++)
        {
            for (integer j = i; j < p; j++)
            {
                resultptr[j + i * n] = tmpptr[idx];
                resultptr[i + j * n] = resultptr[j + i * n];
                idx++;
            }
        }
        for (integer i = 0; i < p; i++)
        {
            for (integer j = p; j < n; j++)
            {
                resultptr[j + i * n] = intretaxptr[idx];
                idx++;
            }
        }
        
        Vector HHR = x.Field("_HHR");
        realdp *HHRptr = HHR.ObtainWritePartialData();
        Vector L(p, p);
        realdp *Lptr = L.ObtainWriteEntireData();
        for(integer i = 0; i < p; i++)
        {
            for(integer j = 0; j < i; j++)
            {
                Lptr[j + i * p] = 0;
            }
            for(integer j = i; j < p; j++)
            {
                Lptr[j + i * p] = HHRptr[i + n * j];
                if(HHRptr[i + n * i] < 0)
                    Lptr[j + i * p] *= -1;
            }
        }
        result->SubmatrixAssignment(0, p - 1, 0, p - 1, L % result->GetSubmatrix(0, p - 1, 0, p - 1));
        
        for (integer i = 0; i < p; i++)
        {
            if(HHRptr[i + n * i] < 0)
                scal_(&p, &GLOBAL::DNONE, resultptr + i, &n);
        }
        
        *result = result->HHRMtp(x.Field("_HHR"), x.Field("_tau"), GLOBAL::N, GLOBAL::L);
        return *result;
	};

	Vector &SymFixedRankQ::EucGradToGradHGZ(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const
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
        Vector XTX(p, p); XTX.AlphaABaddBetaThis(1, x, GLOBAL::T, x, GLOBAL::N, 0);
        *result = egf / XTX;
        return *result;
	};

	Vector &SymFixedRankQ::EucHvToHvHGZ(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const
	{
        Vector egf = x.Field("EGrad");
        
        Vector XTX(p, p); XTX.AlphaABaddBetaThis(1, x, GLOBAL::T, x, GLOBAL::N, 0);
        Vector EGYYinv = egf / XTX;
        Vector YTEGYYinv(p, p); YTEGYYinv.AlphaABaddBetaThis(1, x, GLOBAL::T, EGYYinv, GLOBAL::N, 0);
        Vector EtaTEGYYinv(p, p); EtaTEGYYinv.AlphaABaddBetaThis(1, etax, GLOBAL::T, EGYYinv, GLOBAL::N, 0);
        Vector EtaTY(p, p); EtaTY.AlphaABaddBetaThis(1, etax, GLOBAL::T, x, GLOBAL::N, 0);
        /*xix =  EucHess - 0.5 * Egrad * (Y^T * Y)^{-1} * Y^T * eta - 0.5 * Y * (Y^T * Y)^{-1} Egrad^T * eta + 0.5 * eta * Y^T * Egrad * (Y^T * Y)^{-1}
                 - 0.5 * Y * eta^T * Egrad * (Y^T * Y)^{-1} + 0.5 * eta *(Y^T * Y)^{-1} Egrad^T Y - 0.5 * Egrad * (Y^T * Y)^{-1} * eta^T * Y */
        *result = exix;
        result->AlphaABaddBetaThis(-0.5, EGYYinv, GLOBAL::N, EtaTY, GLOBAL::T, 1);
        result->AlphaABaddBetaThis(-0.5, x, GLOBAL::N, EtaTEGYYinv, GLOBAL::T, 1);
        result->AlphaABaddBetaThis(0.5, etax, GLOBAL::N, YTEGYYinv, GLOBAL::N, 1);
        result->AlphaABaddBetaThis(-0.5, x, GLOBAL::N, EtaTEGYYinv, GLOBAL::N, 1);
        result->AlphaABaddBetaThis(0.5, etax, GLOBAL::N, YTEGYYinv, GLOBAL::T, 1);
        result->AlphaABaddBetaThis(-0.5, EGYYinv, GLOBAL::N, EtaTY, GLOBAL::N, 1); /*result = exix - 0.5 * EGYYinv * EtaTY.GetTranspose() - 0.5 * x * EtaTEGYYinv.GetTranspose() + 0.5 * etax * YTEGYYinv - 0.5 * x * EtaTEGYYinv + 0.5 * etax * YTEGYYinv.GetTranspose() - 0.5 * EGYYinv * EtaTY;*/
        EGYYinv = *result / XTX; /*EGYYinv = result / (x.GetTranspose() * x); EGYYinv is used for temporaray storage*/
        return ExtrProjection(x, EGYYinv, result);
	};

	Vector &SymFixedRankQ::ObtainIntrHGZ(const Variable &x, const Vector &etax, Vector *result)  const
	{
        if(!x.FieldsExist("_HHR"))
        {
            x.HHRDecom();
        }
        
        Vector tmp = etax.HHRMtp(x.Field("_HHR"), x.Field("_tau"), GLOBAL::T, GLOBAL::L);
        Vector HHR = x.Field("_HHR");
        const realdp *HHRptr = HHR.ObtainReadData();
        realdp *tmpptr = tmp.ObtainWritePartialData();
        for (integer i = 0; i < p; i++)
        {
            if(HHRptr[i + n * i] < 0)
                scal_(&p, &GLOBAL::DNONE, tmpptr + i, &n);
        }
        Vector L(p, p);
        realdp *Lptr = L.ObtainWriteEntireData();
        for(integer i = 0; i < p; i++)
        {
            for(integer j = 0; j < i; j++)
            {
                Lptr[j + i * p] = 0;
            }
            for(integer j = i; j < p; j++)
            {
                Lptr[j + i * p] = HHRptr[i + n * j];
                if(HHRptr[i + n * i] < 0)
                    Lptr[j + i * p] *= -1;
            }
        }
        tmp = tmp * L;
        tmpptr = tmp.ObtainWritePartialData(); /*tmp has been updated. we need to reload its pointer*/
        realdp *resultptr = result->ObtainWriteEntireData();
        
        integer idx = 0;
        realdp r2 = static_cast<realdp> (sqrt(2.0));
        
        for (integer i = 0; i < p; i++)
        {
            resultptr[idx] = tmpptr[i + i * n];
            idx++;
        }
        for (integer i = 0; i < p; i++)
        {
            for (integer j = i + 1; j < p; j++)
            {
                resultptr[idx] = r2 * tmpptr[j + i * n];
                idx++;
            }
        }

        for (integer i = 0; i < p; i++)
        {
            for (integer j = p; j < n; j++)
            {
                resultptr[idx] = tmpptr[j + i * n];
                idx++;
            }
        }
        return *result;
	};

	Vector &SymFixedRankQ::ObtainExtrHGZ(const Variable &x, const Vector &intretax, Vector *result) const
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
            resultptr[i + i * n] = intretaxptr[idx];
            idx++;
        }

        for (integer i = 0; i < p; i++)
        {
            for (integer j = i + 1; j < p; j++)
            {
                resultptr[j + i * n] = intretaxptr[idx] / r2;
                resultptr[i + j * n] = intretaxptr[idx] / r2;
                idx++;
            }
        }
        for (integer i = 0; i < p; i++)
        {
            for (integer j = p; j < n; j++)
            {
                resultptr[j + i * n] = intretaxptr[idx];
                idx++;
            }
        }
        for (integer i = 0; i < p; i++)
        {
            if(HHRptr[i + n * i] < 0)
                scal_(&p, &GLOBAL::DNONE, resultptr + i, &n);
        }
        Vector L(p, p);
        realdp *Lptr = L.ObtainWriteEntireData();
        for(integer i = 0; i < p; i++)
        {
            for(integer j = 0; j < i; j++)
            {
                Lptr[j + i * p] = 0;
            }
            for(integer j = i; j < p; j++)
            {
                Lptr[j + i * p] = HHRptr[i + n * j];
                if(HHRptr[i + n * i] < 0)
                    Lptr[j + i * p] *= -1;
            }
        }
        
        *result = *result / L;
        *result = result->HHRMtp(x.Field("_HHR"), x.Field("_tau"), GLOBAL::N, GLOBAL::L);
        return *result;
	};

	Vector &SymFixedRankQ::VectorTransportProj(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const
	{
        if(IsIntrApproach)
        {
            Vector exxix(EMPTYEXTR); ObtainExtr(x, xix, &exxix);
            return ObtainIntr(y, exxix, result);
        }
        
        return ExtrProjection(y, xix, result);
	};

}; /*end of ROPTLIB namespace*/
