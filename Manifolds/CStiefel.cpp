
#include "Manifolds/CStiefel.h"

/*Define the namespace*/
namespace ROPTLIB{

	CStiefel::CStiefel(integer inn, integer inp)
	{
		HasHHR = false;

		metric = CSTIE_EUCLIDEAN;
		retraction = CSTIE_QF;
		VecTran = CSTIE_PARALLELIZATION;
		IsIntrApproach = true;

		n = inn;
		p = inp;
		ExtrinsicDim = 2 * n * p;
		IntrinsicDim = 2 * n * p - p * p;
		name.assign("CStiefel");
		EMPTYEXTR = Vector (n, p, "complex");
		EMPTYINTR = Vector (IntrinsicDim);
	};

	CStiefel::~CStiefel(void)
	{
	};

	void CStiefel::ChooseParamsSet1(void)
	{
		metric = CSTIE_EUCLIDEAN;
		retraction = CSTIE_QF;
		VecTran = CSTIE_PARALLELIZATION;
		IsIntrApproach = true;
		HasHHR = false;
	};

	void CStiefel::ChooseParamsSet2(void)
	{
		metric = CSTIE_EUCLIDEAN;
		retraction = CSTIE_QF;
		VecTran = CSTIE_PROJECTION;
		IsIntrApproach = false;
		HasHHR = false;
	};

	void CStiefel::ChooseParamsSet3(void)
	{
		metric = CSTIE_EUCLIDEAN;
		retraction = CSTIE_POLAR;
		VecTran = CSTIE_PARALLELIZATION;
		IsIntrApproach = true;
		HasHHR = false;
	};

    void CStiefel::ChooseParamsSet4(void)
    {
        metric = CSTIE_EUCLIDEAN;
        retraction = CSTIE_POLAR;
        VecTran = CSTIE_PROJECTION;
        IsIntrApproach = false;
        HasHHR = false;
    };

	void CStiefel::CheckParams(void) const
	{
		std::string StieMetricnames[CSTIEMETRICLENGTH] = { "EUCLIDEAN", "CANONICAL" };
		std::string StieRetractionnames[CSTIERETRACTIONLENGTH] = { "QF", "POLAR", "EXP" };
		std::string StieVectorTransportnames[CSTIEVECTORTRANSPORTLENGTH] = { "PARALLELIZATION", "PARALLELTRANSLATION", "PROJECTION" };
		Manifold::CheckParams();
		printf("%s PARAMETERS:\n", name.c_str());
		printf("n             :%15d,\t", n);
		printf("p             :%15d\n", p);
		printf("metric        :%15s,\t", StieMetricnames[metric].c_str());
		printf("retraction    :%15s\n", StieRetractionnames[retraction].c_str());
		printf("VecTran       :%15s\n", StieVectorTransportnames[VecTran].c_str());
	};

	realdp CStiefel::Metric(const Variable &x, const Vector &etax, const Vector &xix) const
	{
		if (metric == CSTIE_EUCLIDEAN || IsIntrApproach)
			return Manifold::Metric(x, etax, xix);
        else
        if (metric == CSTIE_CANONICAL)
        {
            Vector tmp1(p, p, "complex"); tmp1.AlphaABaddBetaThis(1, etax, GLOBAL::C, x, GLOBAL::N, 0); /* tmp1 = etax.GetTranspose() * x;*/
            Vector tmp2(p, p, "complex"); tmp2.AlphaABaddBetaThis(1, xix, GLOBAL::C, x, GLOBAL::N, 0); /* tmp2 = xix.GetTranspose() * x;*/
            return etax.DotProduct(xix) - 0.5 * tmp1.DotProduct(tmp2);
        }
		printf("Error: Metric has not been done!\n");
		return 0;
	};

    Variable CStiefel::RandominManifold(void) const
    {
        Variable result(n, p, "complex");
        result.RandGaussian();
        result.QRDecom();
        return result.Field("_Q");
    };

	Vector &CStiefel::Projection(const Variable &x, const Vector &etax, Vector *result)const
	{
		if (IsIntrApproach)
			return IntrProjection(x, etax, result);
		else
			return ExtrProjection(x, etax, result);
	};

    Vector &CStiefel::IntrProjection(const Variable &x, const Vector &etax, Vector *result) const
    {
        *result = etax;
        return *result;
    };

    Vector &CStiefel::ExtrProjection(const Variable &x, const Vector &etax, Vector *result) const
    {
        Vector tmp(p, p, "complex"); tmp.AlphaABaddBetaThis(GLOBAL::ZONE, x, GLOBAL::C, etax, GLOBAL::N, GLOBAL::ZZERO); /*tmp = x.GetTranspose() * etax*/
        tmp = (tmp + tmp.GetTranspose()) / 2;
        *result = etax;
        result->AlphaABaddBetaThis(GLOBAL::ZNONE, x, GLOBAL::N, tmp, GLOBAL::N, GLOBAL::ZONE);
        
        return *result;
    };

	Variable &CStiefel::Retraction(const Variable &x, const Vector &etax, Variable *result) const
	{
		if (retraction == CSTIE_QF)
			return qfRetraction(x, etax, result);

		if (retraction == CSTIE_POLAR)
			return PolarRetraction(x, etax, result);

		printf("Error: Retraction has not been done!\n");
        return *result;
	};

    Vector &CStiefel::InvRetraction(const Variable &x, const Variable &y, Vector *result) const
    {
        if (retraction == CSTIE_POLAR)
            return InvPolarRetraction(x, y, result);

        printf("Error: InvRetraction has not been done!\n");
        return Manifold::InvRetraction(x, y, result);
    };

	Vector &CStiefel::coTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
	{
		if (retraction == CSTIE_QF)
			return qfcoTangentVector(x, etax, y, xiy, result);

		if (retraction == CSTIE_POLAR)
			return PolarcoTangentVector(x, etax, y, xiy, result);

		printf("Error: coTangentVector has not been done!\n");
        return Projection(x, xiy, result);
	};

	Vector &CStiefel::DiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir) const
	{
		if (retraction == CSTIE_QF)
			return DiffqfRetraction(x, etax, y, xix, result, IsEtaXiSameDir);

		if (retraction == CSTIE_POLAR)
			return DiffPolarRetraction(x, etax, y, xix, result, IsEtaXiSameDir);

		printf("Error: DiffRetraction has not been done!\n");
        return Manifold::DiffRetraction(x, etax, y, xix, result, IsEtaXiSameDir);
	};

	realdp CStiefel::Beta(const Variable &x, const Vector &etax) const
	{
		if (!HasHHR)
			return 1;

        assert(etax.FieldsExist("beta"));

		/*If the beta has been computed in differentiated retraction, then obtain it.
		Beta should be almost always computed before.*/
        Vector beta = etax.Field("beta");
		const realdp *betav = beta.ObtainReadData();
		return betav[0];
	};

	Vector &CStiefel::VectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const
	{
		if (VecTran == CSTIE_PARALLELIZATION && !HasHHR)
		{
			return Manifold::VectorTransport(x, etax, y, xix, result);
		}

		if (VecTran == CSTIE_PROJECTION && !HasHHR)
			return Projection(y, xix, result);

		if (HasHHR)
			return LCVectorTransport(x, etax, y, xix, result);

		printf("Error: VectorTransport has not been done!\n");
        return Manifold::VectorTransport(x, etax, y, xix, result);
	};

	Vector &CStiefel::InverseVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
	{
		if (VecTran == CSTIE_PARALLELIZATION && !HasHHR)
			return Manifold::InverseVectorTransport(x, etax, y, xiy, result);

		if (VecTran == CSTIE_PROJECTION && !HasHHR)
		{
			printf("CStiefel::InverseVectorTransport: inverse vector transport by projection has not been done!\n");
			return Manifold::InverseVectorTransport(x, etax, y, xiy, result);
		}

		if (HasHHR)
			return LCInverseVectorTransport(x, etax, y, xiy, result);

		printf("Error: CStiefel::InverseVectorTransport has not been done!\n");
        return Manifold::InverseVectorTransport(x, etax, y, xiy, result);
	};

	LinearOPE &CStiefel::HInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, integer start, integer end, LinearOPE *result) const
	{
		if (VecTran == CSTIE_PARALLELIZATION && !HasHHR)
			return Manifold::HInvTran(x, etax, y, Hx, start, end, result);

		if (VecTran == CSTIE_PROJECTION && !HasHHR)
		{
			printf("CStiefel::HInvTran for vector transport by projection has not been done!\n");
            return Manifold::HInvTran(x, etax, y, Hx, start, end, result);
		}

		if (HasHHR)
			return LCHInvTran(x, etax, y, Hx, start, end, result);

		printf("Error: HInvTran has not been done!\n");
        return Manifold::HInvTran(x, etax, y, Hx, start, end, result);
	};

	LinearOPE &CStiefel::TranH(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, integer start, integer end, LinearOPE *result) const
	{
		if (VecTran == CSTIE_PARALLELIZATION && !HasHHR)
			return Manifold::TranH(x, etax, y, Hx, start, end, result);

		if (VecTran == CSTIE_PROJECTION && !HasHHR)
		{
			printf("CStiefel::TranH for vector transport by projection has not been done!\n");
            return Manifold::TranH(x, etax, y, Hx, start, end, result);
		}

		if (HasHHR)
			return LCTranH(x, etax, y, Hx, start, end, result);

		printf("Error: CStiefel::TranH has not been done!\n");
        return Manifold::TranH(x, etax, y, Hx, start, end, result);
	};

	LinearOPE &CStiefel::TranHInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, LinearOPE *result) const
	{
		if (VecTran == CSTIE_PARALLELIZATION && !HasHHR)
			return Manifold::TranHInvTran(x, etax, y, Hx, result);

		if (VecTran == CSTIE_PROJECTION && !HasHHR)
		{
			printf("Warning: CStiefel::TranHInvTran for vector transport by projection has not been done!\n");
            return Manifold::TranHInvTran(x, etax, y, Hx, result);
		}

		if (HasHHR)
			return LCTranHInvTran(x, etax, y, Hx, result);

		printf("Error: CStiefel::TranHInvTran has not been done!\n");
        return Manifold::TranHInvTran(x, etax, y, Hx, result);
	};

	Vector &CStiefel::EucGradToGrad(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const
	{
        if (prob->GetUseHess())
        {
            /*The copy on write is necessary. The reason is that the egf may be from a component in a product of elements.
            Therefore, if CopyOnWrite is not used, then the attached data in x and the product of elements share the same
            memory. This may cause an issue: if the product of elements are released before the attached data in x, then
            release the attached data in x would attempt to delete memory that has been released. This is an error!*/
            Vector EGrad(egf);
            EGrad.CopyOnWrite();
            x.AddToFields("EGrad", EGrad);
        }
		if (metric == CSTIE_EUCLIDEAN)
		{
			return ExtrProjection(x, egf, result);
		}

        /*Canonical metric*/
        Vector tmp(p, p, "complex"); tmp.AlphaABaddBetaThis(1, egf, GLOBAL::C, x, GLOBAL::N, 0);
        *result = egf; result->AlphaABaddBetaThis(-1, x, GLOBAL::N, tmp, GLOBAL::N, 1); /*egf - x * (egf.GetTranspose() * x)*/

        return *result;
	};

	Vector &CStiefel::EucHvToHv(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const
	{
        Vector EGrad = x.Field("EGrad");
		if (metric == CSTIE_EUCLIDEAN)
		{
            Vector tmp(p, p, "complex"); tmp.AlphaABaddBetaThis(1, x, GLOBAL::C, EGrad, GLOBAL::N, 0); /* tmp = x.GetTranspose() * EGrad; */
            tmp = (tmp + tmp.GetTranspose()) / 2;
            *result = exix; result->AlphaABaddBetaThis(-1, etax, GLOBAL::N, tmp, GLOBAL::N, 1);
            return ExtrProjection(x, *result, result);
		}
        /*Canonical metric*/
        Vector tmp(p, p, "complex"); tmp.AlphaABaddBetaThis(1, EGrad, GLOBAL::C, etax, GLOBAL::N, 0); /*tmp = EGrad.GetTranspose() * etax;*/
        Vector tmp2(p, p, "complex"); tmp2.AlphaABaddBetaThis(1, x, GLOBAL::C, EGrad, GLOBAL::N, 0); /*tmp2 = x.GetTranspose() * EGrad;*/

        /*exix - x * (exix.GetTranspose() * x) - x * (tmp - tmp.GetTranspose()) / 2 - (etax * (EGrad.Transpose() * x) - EGrad * (etax.GetTranspose() * x)) / 2
        - 0.5 * ( etax * tmp2 - x * (x.GetTranspose() * etax) * tmp2 )*/
        Vector tmp3(p, p);
        *result = exix;
        tmp3.AlphaABaddBetaThis(1, exix, GLOBAL::C, x, GLOBAL::N, 0);
        result->AlphaABaddBetaThis(-1, x, GLOBAL::N, tmp3, GLOBAL::N, 1); /*exix - x * (exix.GetTranspose() * x)*/
        tmp3 = (tmp - tmp.GetTranspose()) / 2;
        result->AlphaABaddBetaThis(-1, x, GLOBAL::N, tmp3, GLOBAL::N, 1); /*exix - x * (exix.GetTranspose() * x) - x * (tmp - tmp.GetTranspose()) / 2*/
        tmp3.AlphaABaddBetaThis(1, EGrad, GLOBAL::C, x, GLOBAL::N, 0);
        result->AlphaABaddBetaThis(-1, etax, GLOBAL::N, tmp3, GLOBAL::N, 1); /*exix - x * (exix.GetTranspose() * x) - x * (tmp - tmp.GetTranspose()) / 2 - (etax * (EGrad.Transpose() * x)*/
        tmp3.AlphaABaddBetaThis(0.5, etax, GLOBAL::C, x, GLOBAL::N, 0);
        result->AlphaABaddBetaThis(-1, EGrad, GLOBAL::N, tmp3, GLOBAL::N, 1); /*exix - x * (exix.GetTranspose() * x) - x * (tmp - tmp.GetTranspose()) / 2 - (etax * (EGrad.Transpose() * x) - EGrad * (etax.GetTranspose() * x)) / 2*/
        result->AlphaABaddBetaThis(-0.5, etax, GLOBAL::N, tmp2, GLOBAL::N, 1); /*exix - x * (exix.GetTranspose() * x) - x * (tmp - tmp.GetTranspose()) / 2 - (etax * (EGrad.Transpose() * x) - EGrad * (etax.GetTranspose() * x)) / 2 - 0.5 * etax * tmp2*/
        tmp.AlphaABaddBetaThis(1, tmp3, GLOBAL::C, tmp2, GLOBAL::N, 0);
        result->AlphaABaddBetaThis(1, x, GLOBAL::N, tmp, GLOBAL::N, 1); /* exix - x * (exix.GetTranspose() * x) - x * (tmp - tmp.GetTranspose()) / 2 - (etax * (EGrad.Transpose() * x) - EGrad * (etax.GetTranspose() * x)) / 2 - 0.5 * etax * tmp2 + 0.5 * x * (x.GetTranspose() * etax) * tmp2; */

        return *result;
	};

	Vector &CStiefel::ObtainIntr(const Variable &x, const Vector &etax, Vector *result) const
	{
		if (retraction == CSTIE_QF || retraction == CSTIE_POLAR)
			return ObtainIntrHHR(x, etax, result);

        printf("Warning: computing intrinsic representation from extrinsic has not been implemented!\n");
        return ObtainIntrHHR(x, etax, result);
    };

	Vector &CStiefel::ObtainExtr(const Variable &x, const Vector &intretax, Vector *result) const
	{
		if (retraction == CSTIE_QF || retraction == CSTIE_POLAR)
			return ObtainExtrHHR(x, intretax, result);

        printf("Warning: computing extrinsic representation from intrinsic has not been implemented!\n");
        return ObtainExtrHHR(x, intretax, result);
	};

    Vector &CStiefel::ObtainNorVerIntr(const Variable &x, const Vector &etax, Vector *result) const
    {

        Vector xTetax(p, p, "complex"); xTetax.AlphaABaddBetaThis(1, x, GLOBAL::C, etax, GLOBAL::N, 0); /*xTetax  = x.GetTranspose() * etax; */

        realdp *resultptr = result->ObtainWriteEntireData();
        const realdpcomplex *xTetaxptr = (realdpcomplex *) xTetax.ObtainReadData();
        realdp r2 = std::sqrt(2.0);
        integer idx = 0;
        for(integer i = 0; i < p; i++)
        {
            resultptr[idx] = xTetaxptr[i + i * p].r;
            idx++;
        }
        for(integer i = 0; i < p; i++)
        {
            for(integer j = i + 1; j < p; j++)
            {
                resultptr[idx] = (xTetaxptr[i + j * p].r + xTetaxptr[j + i * p].r) / r2;
                idx++;
            }
        }
        for(integer i = 0; i < p; i++)
        {
            for(integer j = i + 1; j < p; j++)
            {
                resultptr[idx] = (xTetaxptr[i + j * p].i - xTetaxptr[j + i * p].i) / r2;
                idx++;
            }
        }
        return *result;
    };

    Vector &CStiefel::ObtainNorVerExtr(const Variable &x, const Vector &intretax, Vector *result) const
    {
        const realdp *intretaxptr = intretax.ObtainReadData();
        Vector xTetax(p, p, "complex");
        realdpcomplex *xTetaxptr = (realdpcomplex *) xTetax.ObtainWriteEntireData();
        realdp r2 = std::sqrt(2.0);
        integer idx = 0;
        for(integer i = 0; i < p; i++)
        {
            xTetaxptr[i + i * p].r = intretaxptr[idx];
            idx++;
        }
        for(integer i = 0; i < p; i++)
        {
            for(integer j = i + 1; j < p; j++)
            {
                xTetaxptr[i + j * p].r = intretaxptr[idx] / r2;
                xTetaxptr[j + i * p].r = xTetaxptr[i + j * p].r;
                idx++;
            }
        }
        for(integer i = 0; i < p; i++)
        {
            for(integer j = i + 1; j < p; j++)
            {
                xTetaxptr[i + j * p].i = intretaxptr[idx] / r2;
                xTetaxptr[j + i * p].i = - xTetaxptr[i + j * p].i;
                idx++;
            }
        }
        *result = x; result->AlphaABaddBetaThis(1, x, GLOBAL::N, xTetax, GLOBAL::N, 0);
        return *result;
    };

	Variable &CStiefel::qfRetraction(const Variable &x, const Vector &etax, Variable *result) const
	{
        Vector exetax(EMPTYEXTR);
        if(IsIntrApproach)
            ObtainExtr(x, etax, &exetax);
        else
            exetax = etax;

        Vector xaddetax = x + exetax;
        xaddetax.QRDecom();
        Vector R = xaddetax.Field("_R");
        *result = xaddetax.Field("_Q");

        realdpcomplex *resultptr = (realdpcomplex *) result->ObtainWritePartialData();
        const realdpcomplex *Rptr = (realdpcomplex *) R.ObtainReadData();
        for(integer i = 0; i < p; i++)
        {
            if(Rptr[i + i * p].r < 0)
            {
                scal_(&n, &GLOBAL::ZNONE, resultptr + i * n, &GLOBAL::IONE);
            }
        }
        result->AddToFields("_HHR", xaddetax.Field("_HHR"));
        result->AddToFields("_tau", xaddetax.Field("_tau"));

        return *result;
	};

    Vector &CStiefel::qfcoTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
    {
        printf("The cotangent vector for the qf retraction has not been implemented!\n");
        return Manifold::coTangentVector(x, etax, y, xiy, result);
    };

    Vector &CStiefel::DiffqfRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir) const
    {
        realdp nxix = std::sqrt(Metric(x, xix, xix));

        Vector exxix(EMPTYEXTR);
        if(IsIntrApproach)
            ObtainExtr(x, xix, &exxix);
        else
            exxix = xix;

        Vector HHR = y.Field("_HHR");
        const realdpcomplex *HHRptr = (realdpcomplex *) HHR.ObtainReadData();
        realdpcomplex *exxixptr = (realdpcomplex *) exxix.ObtainWritePartialData();
        trsm_(GLOBAL::R, GLOBAL::U, GLOBAL::N, GLOBAL::N, &n, &p, &GLOBAL::ZONE, const_cast<realdpcomplex *> (HHRptr), &n, exxixptr, &n);
        for(integer i = 0; i < p; i++)
        {
            if(HHRptr[i + i * n].r < 0)
                scal_(&n, &GLOBAL::ZNONE, exxixptr + i * n, &GLOBAL::IONE);
        }
        Vector YtVRinv(p, p, "complex"); YtVRinv.AlphaABaddBetaThis(1, y, GLOBAL::C, exxix, GLOBAL::N, 0); /*YtVRinv = y.GetTranspose() * exxix; */
        realdpcomplex *YtVRinvptr = (realdpcomplex *) YtVRinv.ObtainWritePartialData();

        for (integer i = 0; i < p; i++)
        {
            YtVRinvptr[i + p * i].r = -YtVRinvptr[i + p * i].r;
            YtVRinvptr[i + p * i].i = 0;
            for (integer j = i + 1; j < p; j++)
            {
                YtVRinvptr[i + p * j].r = -YtVRinvptr[j + p * i].r - YtVRinvptr[i + p * j].r;
                YtVRinvptr[j + p * i].r = 0;
                YtVRinvptr[i + p * j].i = YtVRinvptr[j + p * i].i - YtVRinvptr[i + p * j].i;
                YtVRinvptr[j + p * i].i = 0;
            }
        }
        Vector exresult(exxix); exresult.AlphaABaddBetaThis(1, y, GLOBAL::N, YtVRinv, GLOBAL::N, 1); /*exresult = exxix + y * YtVRinv; */

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

	Variable &CStiefel::PolarRetraction(const Variable &x, const Vector &etax, Variable *result) const
	{
        Vector exetax(EMPTYEXTR);
        if(IsIntrApproach)
            ObtainExtr(x, etax, &exetax);
        else
            exetax = etax;

        exetax.AlphaXaddThis(1, x); /* exetax = x + exetax */
        exetax.SVDDecom();
        *result = exetax.Field("_U") * exetax.Field("_Vt");

        result->AddToFields("_Vt", exetax.Field("_Vt"));
        result->AddToFields("_S", exetax.Field("_S"));
        return *result;
	};

    Vector &CStiefel::InvPolarRetraction(const Variable &x, const Variable &y, Vector *result) const
    {
        Vector A(p, p); A.AlphaABaddBetaThis(1, x, GLOBAL::C, y, GLOBAL::N, 0); /* A = x.GetTranspose() * y; */

        Vector tIp(p, p); tIp.SetToIdentity(); tIp = 2 * tIp;
        Vector soln = tIp.SYL(A, A.GetTranspose());
        *result = x; result->AlphaABaddBetaThis(1, y, GLOBAL::N, soln, GLOBAL::N, -1); /* y * tIp.SYL(A, A.GetTranspose()) - x; */

        return *result;
    };

	Vector &CStiefel::PolarcoTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
	{
		printf("The cotangent vector for the polar retraction has not been implemented!\n");
        return Manifold::coTangentVector(x, etax, y, xiy, result);
	};

	Vector &CStiefel::DiffPolarRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir) const
	{
		if (IsEtaXiSameDir)
		{
            Vector etaxtmp(EMPTYEXTR);
            if(IsIntrApproach)
                ObtainExtr(x, etax, &etaxtmp);
            else
                etaxtmp = etax;
			realdp alpha = sqrt(Metric(x, xix, xix) / Metric(x, etax, etax));
            Vector Vt = y.Field("_Vt");
            Vector S = y.Field("_S");
            Vector tmp(EMPTYEXTR);
            tmp.AlphaABaddBetaThis(1, etaxtmp, GLOBAL::N, Vt, GLOBAL::C, 0); /*tmp = etaxtmp * Vt.GetTranspose();*/
            S = GLOBAL::ZONE / S;
            tmp = S.GetDiagTimesM(tmp, GLOBAL::R) * Vt;
            Vector tmp2(p, p, "complex"); tmp2.AlphaABaddBetaThis(1, tmp, GLOBAL::C, tmp, GLOBAL::N, 0); /*tmp2 = tmp.GetTranspose() * tmp*/
            Vector exresult(tmp); exresult.AlphaABaddBetaThis(-alpha, y, GLOBAL::N, tmp2, GLOBAL::N, alpha); /* exresult = (tmp - y * (tmp.GetTranspose() * tmp)) * alpha; */

            if(IsIntrApproach)
                ObtainIntr(y, exresult, result);
            else
                *result = exresult;

			if (IsEtaXiSameDir && HasHHR)
			{
                realdp nxix = std::sqrt(Metric(x, xix, xix));
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
		}
		printf("Warning: CStiefel::DiffPolarRetraction: The differentiated retraction of the polar retraction has not been implemented!\n");
        return Manifold::DiffRetraction(x, etax, y, xix, result);
	};

	Vector &CStiefel::ObtainIntrHHR(const Variable &x, const Vector &etax, Vector *result) const
	{
        if(!x.FieldsExist("_HHR"))
        {
            x.HHRDecom();
        }
        Vector tmp = etax.HHRMtp(x.Field("_HHR"), x.Field("_tau"), GLOBAL::C, GLOBAL::L);
        Vector HHR = x.Field("_HHR");
        const realdpcomplex *HHRptr = (realdpcomplex *) HHR.ObtainReadData();
        realdpcomplex *tmpptr = (realdpcomplex *) tmp.ObtainWritePartialData();
        for (integer i = 0; i < p; i++)
        {
            if(HHRptr[i + n * i].r < 0)
                scal_(&p, &GLOBAL::ZNONE, tmpptr + i, &n);
        }
        realdp *resultptr = result->ObtainWriteEntireData();

		realdp r2 = static_cast<realdp> (sqrt(2.0));
		integer idx = 0;
        for (integer i = 0; i < p; i++)
        {
            resultptr[idx] = tmpptr[i + i * n].i;
            idx++;
        }
        
		for (integer i = 0; i < p; i++)
		{
			for (integer j = i + 1; j < p; j++)
			{
                resultptr[idx] = r2 * (tmpptr[j + i * n].r - tmpptr[i + j * n].r) / 2;
                idx++;
                resultptr[idx] = r2 * (tmpptr[j + i * n].i + tmpptr[i + j * n].i) / 2;
                idx++;
			}
		}

		for (integer i = 0; i < p; i++)
		{
			for (integer j = p; j < n; j++)
			{
				resultptr[idx] = tmpptr[j + i * n].r;
				idx++;
                resultptr[idx] = tmpptr[j + i * n].i;
                idx++;
			}
		}
        return *result;
	};

	Vector &CStiefel::ObtainExtrHHR(const Variable &x, const Vector &intretax, Vector *result) const
	{
        if(!x.FieldsExist("_HHR"))
        {
            x.HHRDecom();
        }
        realdpcomplex *resultptr = (realdpcomplex *) result->ObtainWriteEntireData();
        const realdp *intretaxptr = intretax.ObtainReadData();
        realdp r2 = static_cast<realdp> (sqrt(2.0));
        integer idx = 0;
        
        for (integer i = 0; i < p; i++)
        {
            resultptr[i + i * n].i = intretaxptr[idx];
            resultptr[i + i * n].r = 0;
            idx++;
        }
        
        for (integer i = 0; i < p; i++)
        {
            for (integer j = i + 1; j < p; j++)
            {
                resultptr[j + i * n].r = intretaxptr[idx] / r2;
                resultptr[i + j * n].r = - resultptr[j + i * n].r;
                idx++;
                resultptr[j + i * n].i = intretaxptr[idx] / r2;
                resultptr[i + j * n].i = resultptr[j + i * n].i;
                idx++;
            }
        }

        for (integer i = 0; i < p; i++)
        {
            for (integer j = p; j < n; j++)
            {
                resultptr[j + i * n].r = intretaxptr[idx];
                idx++;
                resultptr[j + i * n].i = intretaxptr[idx];
                idx++;
            }
        }
        
        Vector HHR = x.Field("_HHR");
        realdpcomplex *HHRptr = (realdpcomplex *) HHR.ObtainWritePartialData();
        for (integer i = 0; i < p; i++)
        {
            if(HHRptr[i + n * i].r < 0)
                scal_(&p, &GLOBAL::ZNONE, resultptr + i, &n);
        }
        (*result) = result->HHRMtp(x.Field("_HHR"), x.Field("_tau"), GLOBAL::N, GLOBAL::L);

        return *result;
	};

	void CStiefel::SetParams(PARAMSMAP params)
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
					break;
				}
			}
		}
	};
}; /*end of ROPTLIB namespace*/
