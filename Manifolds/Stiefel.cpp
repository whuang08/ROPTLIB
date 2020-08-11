
#include "Manifolds/Stiefel.h"

/*Define the namespace*/
namespace ROPTLIB{

	Stiefel::Stiefel(integer inn, integer inp)
	{
		HasHHR = false;

		metric = STIE_EUCLIDEAN;
		retraction = STIE_QF;
		VecTran = STIE_PARALLELIZATION;
		IsIntrApproach = true;

		n = inn;
		p = inp;
		ExtrinsicDim = n * p;
		IntrinsicDim = n * p - p * (p + 1) / 2;
		name.assign("Stiefel");
		EMPTYEXTR = Vector (n, p);
		EMPTYINTR = Vector (IntrinsicDim);
	};

	Stiefel::~Stiefel(void)
	{
	};

	void Stiefel::ChooseParamsSet1(void)
	{
		metric = STIE_EUCLIDEAN;
		retraction = STIE_QF;
		VecTran = STIE_PARALLELIZATION;
		IsIntrApproach = true;
		HasHHR = false;
	};

	void Stiefel::ChooseParamsSet2(void)
	{
		metric = STIE_EUCLIDEAN;
		retraction = STIE_QF;
		VecTran = STIE_PROJECTION;
		IsIntrApproach = false;
		HasHHR = false;
	};

	void Stiefel::ChooseParamsSet3(void)
	{
		metric = STIE_EUCLIDEAN;
		retraction = STIE_POLAR;
		VecTran = STIE_PARALLELIZATION;
		IsIntrApproach = true;
		HasHHR = false;
	};

    void Stiefel::ChooseParamsSet4(void)
    {
        metric = STIE_EUCLIDEAN;
        retraction = STIE_POLAR;
        VecTran = STIE_PROJECTION;
        IsIntrApproach = false;
        HasHHR = false;
    };

    void Stiefel::ChooseParamsSet5(void)
    {
        metric = STIE_EUCLIDEAN;
        retraction = STIE_CAYLEYR;
        VecTran = STIE_CAYLEYVT;
        IsIntrApproach = false;
        HasHHR = false;
    };

	void Stiefel::CheckParams(void) const
	{
		std::string StieMetricnames[STIEMETRICLENGTH] = { "EUCLIDEAN", "CANONICAL" };
		std::string StieRetractionnames[STIERETRACTIONLENGTH] = { "QF", "POLAR", "EXP", "CAYLEYR" };
		std::string StieVectorTransportnames[STIEVECTORTRANSPORTLENGTH] = { "PARALLELIZATION", "RIGGING", "PARALLELTRANSLATION", "PROJECTION", "CAYLEYVT" };
		Manifold::CheckParams();
		printf("%s PARAMETERS:\n", name.c_str());
		printf("n             :%15d,\t", n);
		printf("p             :%15d\n", p);
		printf("metric        :%15s,\t", StieMetricnames[metric].c_str());
		printf("retraction    :%15s\n", StieRetractionnames[retraction].c_str());
		printf("VecTran       :%15s\n", StieVectorTransportnames[VecTran].c_str());
	};

	realdp Stiefel::Metric(const Variable &x, const Vector &etax, const Vector &xix) const
	{
		if (metric == STIE_EUCLIDEAN || IsIntrApproach)
			return Manifold::Metric(x, etax, xix);
        else
        if (metric == STIE_CANONICAL)
        {
            Vector tmp1(p, p); tmp1.AlphaABaddBetaThis(1, etax, GLOBAL::T, x, GLOBAL::N, 0); /* tmp1 = etax.GetTranspose() * x;*/
            Vector tmp2(p, p); tmp2.AlphaABaddBetaThis(1, xix, GLOBAL::T, x, GLOBAL::N, 0); /* tmp2 = xix.GetTranspose() * x;*/
            return etax.DotProduct(xix) - 0.5 * tmp1.DotProduct(tmp2);
        }
		printf("Error: Metric has not been done!\n");
		return 0;
	};

    Variable Stiefel::RandominManifold(void) const
    {
        Variable result(n, p);
        result.RandGaussian();
        result.QRDecom();
        return result.Field("_Q");
    };

	Vector &Stiefel::Projection(const Variable &x, const Vector &etax, Vector *result)const
	{
		if (IsIntrApproach)
			return IntrProjection(x, etax, result);
		else
			return ExtrProjection(x, etax, result);
	};

    Vector &Stiefel::IntrProjection(const Variable &x, const Vector &etax, Vector *result) const
    {
        *result = etax;
        return *result;
    };

    Vector &Stiefel::ExtrProjection(const Variable &x, const Vector &etax, Vector *result) const
    {
        Vector tmp(p, p); tmp.AlphaABaddBetaThis(1, x, GLOBAL::T, etax, GLOBAL::N, 0); /*tmp = x.GetTranspose() * etax*/
        tmp = (tmp + tmp.GetTranspose()) / 2;
        *result = etax;
        result->AlphaABaddBetaThis(-1, x, GLOBAL::N, tmp, GLOBAL::N, 1);
        return *result;
    };

	Variable &Stiefel::Retraction(const Variable &x, const Vector &etax, Variable *result) const
	{
		if (retraction == STIE_QF)
			return qfRetraction(x, etax, result);

		if (retraction == STIE_POLAR)
			return PolarRetraction(x, etax, result);

		if (retraction == STIE_CAYLEYR)
			return CayleyRetraction(x, etax, result);

		printf("Error: Retraction has not been done!\n");
        return *result;
	};

    Vector &Stiefel::InvRetraction(const Variable &x, const Variable &y, Vector *result) const
    {
        if (retraction == STIE_POLAR)
            return InvPolarRetraction(x, y, result);
        
        printf("Error: InvRetraction has not been done!\n");
        return Manifold::InvRetraction(x, y, result);
    };

	Vector &Stiefel::coTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
	{
		if (retraction == STIE_QF)
			return qfcoTangentVector(x, etax, y, xiy, result);

		if (retraction == STIE_POLAR)
			return PolarcoTangentVector(x, etax, y, xiy, result);

		if (retraction == STIE_CAYLEYR)
			return CayleycoTangentVector(x, etax, y, xiy, result);
        
		printf("Error: coTangentVector has not been done!\n");
        return Projection(x, xiy, result);
	};

	Vector &Stiefel::DiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir) const
	{
		if (retraction == STIE_QF)
			return DiffqfRetraction(x, etax, y, xix, result, IsEtaXiSameDir);

		if (retraction == STIE_POLAR)
			return DiffPolarRetraction(x, etax, y, xix, result, IsEtaXiSameDir);

		if (retraction == STIE_CAYLEYR)
			return DiffCayleyRetraction(x, etax, y, xix, result, IsEtaXiSameDir);
        
		printf("Error: DiffRetraction has not been done!\n");
        return Manifold::DiffRetraction(x, etax, y, xix, result, IsEtaXiSameDir);
	};

	realdp Stiefel::Beta(const Variable &x, const Vector &etax) const
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

	Vector &Stiefel::VectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const
	{
		if (VecTran == STIE_PARALLELIZATION && !HasHHR)
		{
			return Manifold::VectorTransport(x, etax, y, xix, result);
		}

		if (VecTran == STIE_PROJECTION && !HasHHR)
			return Projection(y, xix, result);

		if (VecTran == STIE_CAYLEYVT && !HasHHR)
			return CayleyVectorTransport(x, etax, y, xix, result);

		if (HasHHR)
			return LCVectorTransport(x, etax, y, xix, result);
        
		printf("Error: VectorTransport has not been done!\n");
        return Manifold::VectorTransport(x, etax, y, xix, result);
	};

	Vector &Stiefel::InverseVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
	{
		if (VecTran == STIE_PARALLELIZATION && !HasHHR)
			return Manifold::InverseVectorTransport(x, etax, y, xiy, result);

		if (VecTran == STIE_PROJECTION && !HasHHR)
		{
			printf("Stiefel::InverseVectorTransport: inverse vector transport by projection has not been done!\n");
			return Manifold::InverseVectorTransport(x, etax, y, xiy, result);
		}

		if (VecTran == STIE_CAYLEYVT && !HasHHR)
			return CayleyInverseVectorTransport(x, etax, y, xiy, result);

		if (HasHHR)
			return LCInverseVectorTransport(x, etax, y, xiy, result);
        
		printf("Error: Stiefel::InverseVectorTransport has not been done!\n");
        return Manifold::InverseVectorTransport(x, etax, y, xiy, result);
	};

	LinearOPE &Stiefel::HInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, integer start, integer end, LinearOPE *result) const
	{
		if (VecTran == STIE_PARALLELIZATION && !HasHHR)
			return Manifold::HInvTran(x, etax, y, Hx, start, end, result);

		if (VecTran == STIE_PROJECTION && !HasHHR)
		{
			printf("Stiefel::HInvTran for vector transport by projection has not been done!\n");
            return Manifold::HInvTran(x, etax, y, Hx, start, end, result);
		}

		if (VecTran == STIE_CAYLEYVT && !HasHHR)
		{
			printf("Stiefel::HInvTran for Cayley vector transport has not been done!\n");
            return Manifold::HInvTran(x, etax, y, Hx, start, end, result);
		}

		if (HasHHR)
			return LCHInvTran(x, etax, y, Hx, start, end, result);

		printf("Error: HInvTran has not been done!\n");
        return Manifold::HInvTran(x, etax, y, Hx, start, end, result);
	};

	LinearOPE &Stiefel::TranH(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, integer start, integer end, LinearOPE *result) const
	{
		if (VecTran == STIE_PARALLELIZATION && !HasHHR)
			return Manifold::TranH(x, etax, y, Hx, start, end, result);

		if (VecTran == STIE_PROJECTION && !HasHHR)
		{
			printf("Stiefel::TranH for vector transport by projection has not been done!\n");
            return Manifold::TranH(x, etax, y, Hx, start, end, result);
		}

		if (VecTran == STIE_CAYLEYVT && !HasHHR)
		{
			printf("Stiefel::TranH for Cayley vector transport has not been done!\n");
            return Manifold::TranH(x, etax, y, Hx, start, end, result);
		}

		if (HasHHR)
			return LCTranH(x, etax, y, Hx, start, end, result);

		printf("Error: Stiefel::TranH has not been done!\n");
        return Manifold::TranH(x, etax, y, Hx, start, end, result);
	};

	LinearOPE &Stiefel::TranHInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, LinearOPE *result) const
	{
		if (VecTran == STIE_PARALLELIZATION && !HasHHR)
			return Manifold::TranHInvTran(x, etax, y, Hx, result);

		if (VecTran == STIE_PROJECTION && !HasHHR)
		{
			printf("Warning: Stiefel::TranHInvTran for vector transport by projection has not been done!\n");
            return Manifold::TranHInvTran(x, etax, y, Hx, result);
		}

		if (VecTran == STIE_CAYLEYVT && !HasHHR)
		{
			printf("Warning: Stiefel::TranHInvTran for Cayley vector transport has not been done!\n");
            return Manifold::TranHInvTran(x, etax, y, Hx, result);
		}

		if (HasHHR)
			return LCTranHInvTran(x, etax, y, Hx, result);

		printf("Error: Stiefel::TranHInvTran has not been done!\n");
        return Manifold::TranHInvTran(x, etax, y, Hx, result);
	};

	Vector &Stiefel::EucGradToGrad(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const
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
		if (metric == STIE_EUCLIDEAN)
		{
			return ExtrProjection(x, egf, result);
		}

        /*Canonical metric*/
        Vector tmp(p, p); tmp.AlphaABaddBetaThis(1, egf, GLOBAL::T, x, GLOBAL::N, 0);
        *result = egf; result->AlphaABaddBetaThis(-1, x, GLOBAL::N, tmp, GLOBAL::N, 1); /*egf - x * (egf.GetTranspose() * x)*/
        
        return *result;
	};

	Vector &Stiefel::EucHvToHv(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const
	{
        Vector EGrad = x.Field("EGrad");
		if (metric == STIE_EUCLIDEAN)
		{
            Vector tmp(p, p); tmp.AlphaABaddBetaThis(1, x, GLOBAL::T, EGrad, GLOBAL::N, 0); /* tmp = x.GetTranspose() * EGrad; */
            tmp = (tmp + tmp.GetTranspose()) / 2;
            *result = exix; result->AlphaABaddBetaThis(-1, etax, GLOBAL::N, tmp, GLOBAL::N, 1);
            return ExtrProjection(x, *result, result);
		}
        /*Canonical metric*/
        Vector tmp(p, p); tmp.AlphaABaddBetaThis(1, EGrad, GLOBAL::T, etax, GLOBAL::N, 0); /*tmp = EGrad.GetTranspose() * etax;*/
        Vector tmp2(p, p); tmp2.AlphaABaddBetaThis(1, x, GLOBAL::T, EGrad, GLOBAL::N, 0); /*tmp2 = x.GetTranspose() * EGrad;*/
        
        /*exix - x * (exix.GetTranspose() * x) - x * (tmp - tmp.GetTranspose()) / 2 - (etax * (EGrad.Transpose() * x) - EGrad * (etax.GetTranspose() * x)) / 2
        - 0.5 * ( etax * tmp2 - x * (x.GetTranspose() * etax) * tmp2 )*/
        Vector tmp3(p, p);
        *result = exix;
        tmp3.AlphaABaddBetaThis(1, exix, GLOBAL::T, x, GLOBAL::N, 0);
        result->AlphaABaddBetaThis(-1, x, GLOBAL::N, tmp3, GLOBAL::N, 1); /*exix - x * (exix.GetTranspose() * x)*/
        tmp3 = (tmp - tmp.GetTranspose()) / 2;
        result->AlphaABaddBetaThis(-1, x, GLOBAL::N, tmp3, GLOBAL::N, 1); /*exix - x * (exix.GetTranspose() * x) - x * (tmp - tmp.GetTranspose()) / 2*/
        tmp3.AlphaABaddBetaThis(1, EGrad, GLOBAL::T, x, GLOBAL::N, 0);
        result->AlphaABaddBetaThis(-1, etax, GLOBAL::N, tmp3, GLOBAL::N, 1); /*exix - x * (exix.GetTranspose() * x) - x * (tmp - tmp.GetTranspose()) / 2 - (etax * (EGrad.Transpose() * x)*/
        tmp3.AlphaABaddBetaThis(0.5, etax, GLOBAL::T, x, GLOBAL::N, 0);
        result->AlphaABaddBetaThis(-1, EGrad, GLOBAL::N, tmp3, GLOBAL::N, 1); /*exix - x * (exix.GetTranspose() * x) - x * (tmp - tmp.GetTranspose()) / 2 - (etax * (EGrad.Transpose() * x) - EGrad * (etax.GetTranspose() * x)) / 2*/
        result->AlphaABaddBetaThis(-0.5, etax, GLOBAL::N, tmp2, GLOBAL::N, 1); /*exix - x * (exix.GetTranspose() * x) - x * (tmp - tmp.GetTranspose()) / 2 - (etax * (EGrad.Transpose() * x) - EGrad * (etax.GetTranspose() * x)) / 2 - 0.5 * etax * tmp2*/
        tmp.AlphaABaddBetaThis(1, tmp3, GLOBAL::T, tmp2, GLOBAL::N, 0);
        result->AlphaABaddBetaThis(1, x, GLOBAL::N, tmp, GLOBAL::N, 1); /* exix - x * (exix.GetTranspose() * x) - x * (tmp - tmp.GetTranspose()) / 2 - (etax * (EGrad.Transpose() * x) - EGrad * (etax.GetTranspose() * x)) / 2 - 0.5 * etax * tmp2 + 0.5 * x * (x.GetTranspose() * etax) * tmp2; */
        
        return *result;
	};

	Vector &Stiefel::ObtainIntr(const Variable &x, const Vector &etax, Vector *result) const
	{
		if (retraction == STIE_QF || retraction == STIE_POLAR)
			return ObtainIntrHHR(x, etax, result);
        
        printf("Warning: computing intrinsic representation from extrinsic has not been implemented!\n");
        return ObtainIntrHHR(x, etax, result);
    };

	Vector &Stiefel::ObtainExtr(const Variable &x, const Vector &intretax, Vector *result) const
	{
		if (retraction == STIE_QF || retraction == STIE_POLAR)
			return ObtainExtrHHR(x, intretax, result);
		
        printf("Warning: computing extrinsic representation from intrinsic has not been implemented!\n");
        return ObtainExtrHHR(x, intretax, result);
	};

    Vector &Stiefel::ObtainNorVerIntr(const Variable &x, const Vector &etax, Vector *result) const
    {
        
        Vector xTetax(p, p); xTetax.AlphaABaddBetaThis(1, x, GLOBAL::T, etax, GLOBAL::N, 0); /*xTetax  = x.GetTranspose() * etax; */
        
        realdp *resultptr = result->ObtainWriteEntireData();
        const realdp *xTetaxptr = xTetax.ObtainReadData();
        realdp r2 = std::sqrt(2.0);
        integer idx = 0;
        for(integer i = 0; i < p; i++)
        {
            resultptr[idx] = xTetaxptr[i + i * p];
            idx++;
        }
        for(integer i = 0; i < p; i++)
        {
            for(integer j = i + 1; j < p; j++)
            {
                resultptr[idx] = (xTetaxptr[i + j * p] + xTetaxptr[j + i * p]) / r2;
                idx++;
            }
        }
        return *result;
    };

    Vector &Stiefel::ObtainNorVerExtr(const Variable &x, const Vector &intretax, Vector *result) const
    {
        const realdp *intretaxptr = intretax.ObtainReadData();
        Vector xTetax(p, p);
        realdp *xTetaxptr = xTetax.ObtainWriteEntireData();
        realdp r2 = std::sqrt(2.0);
        integer idx = 0;
        for(integer i = 0; i < p; i++)
        {
            xTetaxptr[i + i * p] = intretaxptr[idx];
            idx++;
        }
        for(integer i = 0; i < p; i++)
        {
            for(integer j = i + 1; j < p; j++)
            {
                xTetaxptr[i + j * p] = intretaxptr[idx] / r2;
                xTetaxptr[j + i * p] = xTetaxptr[i + j * p];
                idx++;
            }
        }
        *result = x; result->AlphaABaddBetaThis(1, x, GLOBAL::N, xTetax, GLOBAL::N, 0);
        return *result;
    };

	Variable &Stiefel::qfRetraction(const Variable &x, const Vector &etax, Variable *result) const
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
        
        realdp *resultptr = result->ObtainWritePartialData();
        const realdp *Rptr = R.ObtainReadData();
        for(integer i = 0; i < p; i++)
        {
            if(Rptr[i + i * p] < 0)
            {
                scal_(&n, &GLOBAL::DNONE, resultptr + i * n, &GLOBAL::IONE);
            }
        }
        result->AddToFields("_HHR", xaddetax.Field("_HHR"));
        result->AddToFields("_tau", xaddetax.Field("_tau"));
        
        return *result;
	};

    Vector &Stiefel::qfcoTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
    {
        Vector exxiy(EMPTYEXTR);
        if(IsIntrApproach)
            ObtainExtr(y, xiy, &exxiy);
        else
            exxiy = xiy;
        
        Vector ytxiy(p, p); ytxiy.AlphaABaddBetaThis(1, y, GLOBAL::T, exxiy, GLOBAL::N, 0); /* ytxiy = y.GetTranspose() * exxiy;*/
        realdp *ytxiyptr = ytxiy.ObtainWritePartialData();
        for (integer i = 0; i < p; i++)
        {
            for (integer j = i; j < p; j++)
            {
                ytxiyptr[i + j * p] = -ytxiyptr[i + j * p];
            }
        }
        
        Vector exresult = exxiy; exresult.AlphaABaddBetaThis(1, y, GLOBAL::N, ytxiy, GLOBAL::N, 1); /*exresult = y * ytxiy + exxiy;*/
        
        realdp *exresultptr = exresult.ObtainWritePartialData();
        Vector HHR = y.Field("_HHR");
        const realdp *HHRptr = HHR.ObtainReadData();
        for (integer i = 0; i < p; i++)
        {
            if(HHRptr[i + i * n] < 0)
                scal_(&n, &GLOBAL::DNONE, exresultptr + i * n, &GLOBAL::IONE);
        }
        
        trsm_(GLOBAL::R, GLOBAL::U, GLOBAL::T, GLOBAL::N, &n, &p, &GLOBAL::DONE, const_cast<realdp *> (HHRptr), &n, exresultptr, &n);
        ExtrProjection(x, exresult, &exresult);
        
        if (IsIntrApproach)
        {
            return ObtainIntr(x, exresult, result);
        }
        *result = exresult;
        return *result;
    };

    Vector &Stiefel::DiffqfRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir) const
    {
        realdp nxix = std::sqrt(Metric(x, xix, xix));

        Vector exxix(EMPTYEXTR);
        if(IsIntrApproach)
            ObtainExtr(x, xix, &exxix);
        else
            exxix = xix;
        
        Vector HHR = y.Field("_HHR");
        const realdp *HHRptr = HHR.ObtainReadData();
        realdp *exxixptr = exxix.ObtainWritePartialData();
        trsm_(GLOBAL::R, GLOBAL::U, GLOBAL::N, GLOBAL::N, &n, &p, &GLOBAL::DONE, const_cast<realdp *> (HHRptr), &n, exxixptr, &n);
        for(integer i = 0; i < p; i++)
        {
            if(HHRptr[i + i * n] < 0)
                scal_(&n, &GLOBAL::DNONE, exxixptr + i * n, &GLOBAL::IONE);
        }
        Vector YtVRinv(p, p); YtVRinv.AlphaABaddBetaThis(1, y, GLOBAL::T, exxix, GLOBAL::N, 0); /*YtVRinv = y.GetTranspose() * exxix; */
        realdp *YtVRinvptr = YtVRinv.ObtainWritePartialData();
        
        for (integer i = 0; i < p; i++)
        {
            YtVRinvptr[i + p * i] = -YtVRinvptr[i + p * i];
            for (integer j = i + 1; j < p; j++)
            {
                YtVRinvptr[i + p * j] = -YtVRinvptr[j + p * i] - YtVRinvptr[i + p * j];
                YtVRinvptr[j + p * i] = 0;
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

	Variable &Stiefel::PolarRetraction(const Variable &x, const Vector &etax, Variable *result) const
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

    Vector &Stiefel::InvPolarRetraction(const Variable &x, const Variable &y, Vector *result) const
    {
        Vector A(p, p); A.AlphaABaddBetaThis(1, x, GLOBAL::T, y, GLOBAL::N, 0); /* A = x.GetTranspose() * y; */
        
        Vector tIp(p, p); tIp.SetToIdentity(); tIp = 2 * tIp;
        Vector soln = tIp.SYL(A, A.GetTranspose());
        *result = x; result->AlphaABaddBetaThis(1, y, GLOBAL::N, soln, GLOBAL::N, -1); /* y * tIp.SYL(A, A.GetTranspose()) - x; */
        
        return *result;
    };

	Vector &Stiefel::PolarcoTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
	{
		printf("The cotangent vector for the polar retraction has not been implemented!\n");
        return Manifold::coTangentVector(x, etax, y, xiy, result);
	};

	Vector &Stiefel::DiffPolarRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir) const
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
            tmp.AlphaABaddBetaThis(1, etaxtmp, GLOBAL::N, Vt, GLOBAL::T, 0); /*tmp = etaxtmp * Vt.GetTranspose();*/
            S = 1 / S;
            tmp = S.GetDiagTimesM(tmp, GLOBAL::R) * Vt;
            Vector tmp2(p, p); tmp2.AlphaABaddBetaThis(1, tmp, GLOBAL::T, tmp, GLOBAL::N, 0); /*tmp2 = tmp.GetTranspose() * tmp*/
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
		printf("Warning: Stiefel::DiffPolarRetraction: The differentiated retraction of the polar retraction has not been implemented!\n");
        return Manifold::DiffRetraction(x, etax, y, xix, result);
	};

	Vector &Stiefel::ObtainIntrHHR(const Variable &x, const Vector &etax, Vector *result) const
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
        
        realdp *resultptr = result->ObtainWriteEntireData();
        
		realdp r2 = static_cast<realdp> (sqrt(2.0));
		integer idx = 0;
		for (integer i = 0; i < p; i++)
		{
			for (integer j = i + 1; j < p; j++)
			{
				resultptr[idx] = r2 * (tmpptr[j + i * n] - tmpptr[i + j * n]) / 2;
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

	Vector &Stiefel::ObtainExtrHHR(const Variable &x, const Vector &intretax, Vector *result) const
	{
        if(!x.FieldsExist("_HHR"))
        {
            x.HHRDecom();
        }
        realdp *resultptr = result->ObtainWriteEntireData();
        const realdp *intretaxptr = intretax.ObtainReadData();
        realdp r2 = static_cast<realdp> (sqrt(2.0));
        integer idx = 0;
        for (integer i = 0; i < p; i++)
        {
            resultptr[i + i * n] = 0;
            for (integer j = i + 1; j < p; j++)
            {
                resultptr[j + i * n] = intretaxptr[idx] / r2;
                resultptr[i + j * n] = -resultptr[j + i * n];
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
        for (integer i = 0; i < p; i++)
        {
            if(HHRptr[i + n * i] < 0)
                scal_(&p, &GLOBAL::DNONE, resultptr + i, &n);
        }
        (*result) = result->HHRMtp(x.Field("_HHR"), x.Field("_tau"), GLOBAL::N, GLOBAL::L);
        
        return *result;
	};

	Variable &Stiefel::CayleyRetraction(const Variable &x, const Vector &etax, Variable *result) const
	{/* assume extrinsic representation is used for etax */
        Vector exetax(EMPTYEXTR);
        if(IsIntrApproach)
            ObtainExtr(x, etax, &exetax);
        else
            exetax = etax;
        
        Vector Mk2(2 * p, 2 * p), Mk3(2 * p, p), LUP(4 * p * p + 2 * p), U(n, 2 * p);
        realdp *Mk2ptr = Mk2.ObtainWriteEntireData();
        realdp *Mk3ptr = Mk3.ObtainWriteEntireData();
        realdp *LUPptr = LUP.ObtainWriteEntireData();
        realdp *Uptr = U.ObtainWriteEntireData();
        Vector Mk211(p, p); Mk211.AlphaABaddBetaThis(0.5, x, GLOBAL::T, etax, GLOBAL::N, 0); /* Mk211 = 0.5 * x.GetTranspose() * etax; */
        
		const realdp *xptr = x.ObtainReadData();
		const realdp *etaxptr = etax.ObtainReadData();

		integer P = p, N = n, P2 = p * 2;
		realdp half = 0.5;
		/*1-1 block of Mk2*/
		gemm_(GLOBAL::T, GLOBAL::N, &P, &P, &N, &half, const_cast<realdp *> (xptr), &N, const_cast<realdp *> (etaxptr), &N, &GLOBAL::DZERO, Mk2ptr, &P2);
		/*1-2 block of Mk2*/
		for (integer i = 0; i < p; i++)
		{
			for (integer j = p; j < P2; j++)
			{
				Mk2ptr[i + j * P2] = 0;
			}
			Mk2ptr[i + (i + p) * P2] = 1;
		}
		/*2-1 block of Mk2*/
		realdp three = 3.0;
		gemm_(GLOBAL::T, GLOBAL::N, &P, &P, &P, &three, Mk2ptr, &P2, Mk2ptr, &P2, &GLOBAL::DZERO, Mk2ptr + p, &P2);
		gemm_(GLOBAL::T, GLOBAL::N, &P, &P, &N, &GLOBAL::DNONE, const_cast<realdp *> (etaxptr), &N, const_cast<realdp *> (etaxptr), &N, &GLOBAL::DONE, Mk2ptr + p, &P2);
		/*2-2 block of Mk2*/
		for (integer i = 0; i < p; i++)
		{
			for (integer j = 0; j < p; j++)
			{
				Mk2ptr[i + p + (j + p) * P2] = -Mk2ptr[j + i * P2];
			}
		}
		realdp *Mk1ptr = Mk2ptr + P2 * p;
		/*for computing LUP*/
		integer Psquare4 = p * p * 4;
		realdp nhalf = -0.5;
		for (integer i = 0; i < 4 * p * p; i++)
			LUPptr[i] = 0;
		axpy_(&Psquare4, &nhalf, Mk2ptr, &GLOBAL::IONE, LUPptr, &GLOBAL::IONE);
		for (integer i = 0; i < P2; i++)
			LUPptr[i + i * P2] += 1;
		integer info;
		integer *Perm = new integer[P2];
		/* LU decomposion for LUPpter, LUPpter = P * L * U, L and U are stored in LUPpter, the permutation matrix is in Perm
		 details: www.netlib.org/lapack/explore-html/d3/d6a/getrf_8f.html */
		getrf_(&P2, &P2, LUPptr, &P2, Perm, &info);
		for (integer i = 0; i < P2; i++)
			LUPptr[4 * p * p + i] = static_cast<realdp> (Perm[i]);
		/*Compute Mk3*/
		integer length = 2 * p * p;
		copy_(&length, Mk1ptr, &GLOBAL::IONE, Mk3ptr, &GLOBAL::IONE);
		/*Solve the linear system*/
		getrs_(GLOBAL::N, &P2, &P, LUPptr, &P2, Perm, Mk3ptr, &P2, &info);
		if (info != 0)
			printf("Warning: getrs in Stiefel::CayleyRetraction failed!\n");
		delete[] Perm;
		/*Compute U*/
		length = n * p;
		copy_(&length, const_cast<realdp *>(etaxptr), &GLOBAL::IONE, Uptr, &GLOBAL::IONE);
		copy_(&length, const_cast<realdp *>(xptr), &GLOBAL::IONE, Uptr + n * p, &GLOBAL::IONE);
		gemm_(GLOBAL::N, GLOBAL::N, &N, &P, &P, &GLOBAL::DNONE, const_cast<realdp *> (xptr), &N, Mk2ptr, &P2, &GLOBAL::DONE, Uptr, &N);
		/*compute the result: x + U * Mk3*/

        *result = x;
		realdp *resultptr = result->ObtainWritePartialData();
		gemm_(GLOBAL::N, GLOBAL::N, &N, &P, &P2, &GLOBAL::DONE, Uptr, &N, Mk3ptr, &P2, &GLOBAL::DONE, resultptr, &N);

        x.AddToFields("Mk2", Mk2);
        x.AddToFields("Mk3", Mk3);
        x.AddToFields("LUP", LUP);
        x.AddToFields("U", U);
        return *result;
	};

	Vector &Stiefel::CayleycoTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
	{
		printf("The cotangent vector for the Cayley retraction has not been implemented!\n");
        return Manifold::coTangentVector(x, etax, y, xiy, result);
	};

	Vector &Stiefel::DiffCayleyRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir) const
	{
		if (IsEtaXiSameDir)
		{
			const realdp *Mk2ptr = x.Field("Mk2").ObtainReadData();
			const realdp *Mk3ptr = x.Field("Mk3").ObtainReadData();
			const realdp *LUPptr = x.Field("LUP").ObtainReadData();
			const realdp *Uptr = x.Field("U").ObtainReadData();

			realdp *Mk2Mk3 = new realdp[4 * p * p];
			realdp *tmp = Mk2Mk3 + 2 * p * p;
			integer P2 = p * 2, P = p, N = n, length = 2 * p * p;
			/*compute Mk2 * Mk3 */
			gemm_(GLOBAL::N, GLOBAL::N, &P2, &P, &P2, &GLOBAL::DONE, const_cast<realdp *> (Mk2ptr), &P2, const_cast<realdp *> (Mk3ptr), &P2, &GLOBAL::DZERO, Mk2Mk3, &P2);
			/*tmp <-- Mk1*/
			copy_(&length, const_cast<realdp *> (Mk2ptr + 2 * p * p), &GLOBAL::IONE, tmp, &GLOBAL::IONE);
			realdp half = 0.5;
			/*tmp <-- Mk1 + 0.5 * Mk2 * Mk3 */
			axpy_(&length, &half, Mk2Mk3, &GLOBAL::IONE, tmp, &GLOBAL::IONE);

			/*solve (1 - 0.5Mk2)^{-1} Mk2Mk3*/
			integer info, *perm = new integer[2 * p];
			for (integer i = 0; i < 2 * p; i++)
				perm[i] = static_cast<integer> (LUPptr[4 * p * p + i]);

			/*Solve the linear system*/
			getrs_(GLOBAL::N, &P2, &P, const_cast<realdp *>(LUPptr), &P2, perm, Mk2Mk3, &P2, &info);
			if (info != 0)
				printf("Warning: getrs in Stiefel::DiffCayleyRetraction failed!\n");
			delete[] perm;
			/*tmp <--Mk1 + 0.5 * Mk2 * Mk3 + 0.5 (1 - 0.5Mk2)^{-1} Mk2Mk3 */
			axpy_(&length, &half, Mk2Mk3, &GLOBAL::IONE, tmp, &GLOBAL::IONE);

            Vector exresult(n, p);
			realdp *exresultptr = exresult.ObtainWriteEntireData();
			gemm_(GLOBAL::N, GLOBAL::N, &N, &P, &P2, &GLOBAL::DONE, const_cast<realdp *> (Uptr), &N, tmp, &P2, &GLOBAL::DZERO, exresultptr, &N);

			delete[] Mk2Mk3;
            if(IsIntrApproach)
                ObtainIntr(x, exresult, result);
            else
                *result = exresult;
            
			ScalarTimesVector(x, sqrt(Metric(x, xix, xix) / Metric(x, etax, etax)), *result, result);

            if (IsEtaXiSameDir && HasHHR)
            {
                Vector beta(3);
                realdp *betaptr = beta.ObtainWriteEntireData();
                realdp EtatoXi = sqrt(Metric(x, etax, etax) / Metric(x, xix, xix));
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
        
		printf("Warning: Stiefel::DiffCayleyRetraction: The differentiated retraction of the Cayley retraction has not been implemented!\n");
        return Manifold::DiffRetraction(x, etax, y, xix, result);
	};

	Vector &Stiefel::CayleyVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const
	{
        Vector exxix(EMPTYEXTR);
        if(IsIntrApproach)
            ObtainExtr(x, xix, &exxix);
        else
            exxix = xix;
        
        const realdp *LUPptr = x.Field("LUP").ObtainReadData();
        const realdp *Uptr = x.Field("U").ObtainReadData();
        
		realdp *V = new realdp[n * 2 * p + 2 * p * p];
		realdp *Vtxix = V + n * 2 * p;
		integer length = n * p;
		copy_(&length, const_cast<realdp *> (Uptr + n * p), &GLOBAL::IONE, V, &GLOBAL::IONE);
		copy_(&length, const_cast<realdp *> (Uptr), &GLOBAL::IONE, V + n * p, &GLOBAL::IONE);
		scal_(&length, &GLOBAL::DNONE, V + n * p, &GLOBAL::IONE);
		const realdp *xixptr = xix.ObtainReadData();
		integer N = n, P = p, P2 = p * 2;
		gemm_(GLOBAL::T, GLOBAL::N, &P2, &P, &N, &GLOBAL::DONE, V, &N, const_cast<realdp *> (xixptr), &N, &GLOBAL::DZERO, Vtxix, &P2);

		integer info, *Perm = new integer[P2];
		for (integer i = 0; i < P2; i++)
			Perm[i] = static_cast<integer> (LUPptr[4 * p * p + i]);
		/*Solve the linear system*/
		getrs_(GLOBAL::N, &P2, &P, const_cast<realdp *>(LUPptr), &P2, Perm, Vtxix, &P2, &info);
		if (info != 0)
			printf("Warning: getrs in Stiefel::DiffCayleyRetraction failed!\n");
		delete[] Perm;

        Vector exresult(xix);
		realdp *exresultptr = exresult.ObtainWritePartialData();
		gemm_(GLOBAL::N, GLOBAL::N, &N, &P, &P2, &GLOBAL::DONE, const_cast<realdp *> (Uptr), &N, Vtxix, &P2, &GLOBAL::DONE, exresultptr, &N);
		delete[] V;
        if(IsIntrApproach)
            return ObtainIntr(y, exresult, result);
        
        *result = exresult;
        return *result;
	};

	Vector &Stiefel::CayleyInverseVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
	{
        Vector exxiy(EMPTYEXTR);
        if(IsIntrApproach)
            ObtainExtr(y, xiy, &exxiy);
        else
            exxiy = xiy;
        
        const realdp *LUPptr = x.Field("LUP").ObtainReadData();
        const realdp *Uptr = x.Field("U").ObtainReadData();
        
		realdp *V = new realdp[n * 2 * p + 2 * p * p];
		realdp *Utxiy = V + n * 2 * p;
		integer length = n * p;
		copy_(&length, const_cast<realdp *> (Uptr + n * p), &GLOBAL::IONE, V, &GLOBAL::IONE);
		copy_(&length, const_cast<realdp *> (Uptr), &GLOBAL::IONE, V + n * p, &GLOBAL::IONE);
		scal_(&length, &GLOBAL::DNONE, V + n * p, &GLOBAL::IONE);
		const realdp *xiyptr = xiy.ObtainReadData();
		integer N = n, P = p, P2 = p * 2;
		gemm_(GLOBAL::T, GLOBAL::N, &P2, &P, &N, &GLOBAL::DONE, const_cast<realdp *> (Uptr), &N, const_cast<realdp *> (xiyptr), &N, &GLOBAL::DZERO, Utxiy, &P2);

		integer info, *Perm = new integer[P2];
		for (integer i = 0; i < P2; i++)
			Perm[i] = static_cast<integer> (LUPptr[4 * p * p + i]);
		/*Solve the linear system*/
		getrs_(GLOBAL::T, &P2, &P, const_cast<realdp *>(LUPptr), &P2, Perm, Utxiy, &P2, &info);
		if (info != 0)
			printf("Warning: getrs in Stiefel::DiffCayleyRetraction failed!\n");
		delete[] Perm;

        Vector exresult(xiy);
		realdp *exresultptr = exresult.ObtainWritePartialData();
		gemm_(GLOBAL::N, GLOBAL::N, &N, &P, &P2, &GLOBAL::DONE, V, &N, Utxiy, &P2, &GLOBAL::DONE, exresultptr, &N);
		delete[] V;
        if(IsIntrApproach)
            return ObtainIntr(x, exresult, result);
        
        *result = exresult;
        return *result;
	};

	void Stiefel::SetParams(PARAMSMAP params)
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
				default:
					break;
				}
			}
		}
	};
}; /*end of ROPTLIB namespace*/
