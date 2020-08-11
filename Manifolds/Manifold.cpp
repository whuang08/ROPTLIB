
#include "Manifolds/Manifold.h"

/*Define the namespace*/
namespace ROPTLIB{

	Manifold::~Manifold(void)
	{
	};

	realdp Manifold::Metric(const Variable &x, const Vector &etax, const Vector &xix) const
	{
    #ifdef CHECKMANIFOLDOVERRIDDEN
            printf("Manifold::Metric has not been overridden!\n");
    #endif
        return etax.DotProduct(xix);
	};

	Vector &Manifold::LinearOPEEta(const Variable &x, const LinearOPE &Hx, const Vector &etax, Vector *result) const
	{
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::LinearOPEEta has not been overridden!\n");
        #endif

        /*even if etax is a complex element, the real reparameterizations are used. In other words, we viewed it is
        2 n dimension real data. We have to convert it to real element to apply the linear opeartor. */
        bool etaxiscomplex = etax.Getiscomplex();
        etax.Setiscomplex(false);
        Vector etaxreshape = etax; etaxreshape.Reshape(etax.Getlength());
        result->AlphaABaddBetaThis(1, Hx, GLOBAL::N, etaxreshape, GLOBAL::N, 0);
        etax.Setiscomplex(etaxiscomplex);
        result->Setiscomplex(etaxiscomplex);
        result->Reshape(etax.Getrow(), etax.Getcol(), etax.Getnum());
        return *result;
	};

    Vector &Manifold::Projection(const Variable &x, const Vector &etax, Vector *result) const
    {
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::Projection has not been overridden!\n");
        #endif
        if (!IsIntrApproach)
        {
            return this->ExtrProjection(x, etax, result);
        }
        else
        {
            return this->IntrProjection(x, etax, result);
        }
    };

    Vector &Manifold::IntrProjection(const Variable &x, const Vector &etax, Vector *result) const
    {
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::IntrProjection has not been overridden!\n");
        #endif
        *result = etax;
        return *result;
    };

    Vector &Manifold::ExtrProjection(const Variable &x, const Vector &etax, Vector *result) const
    {
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::ExtrProjection has not been overridden!\n");
        #endif
        *result = etax;
        return *result;
    };

	Vector &Manifold::ScalarTimesVector(const Variable &x, const realdp &scalar, const Vector &etax, Vector *result) const
	{
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::ScalarTimesVector has not been overridden!\n");
        #endif
        *result = etax;
        result->ScalarTimesThis(scalar);
        return *result;
	};

    Vector &Manifold::ScalarVectorAddVector(const Variable &x, const realdp &scalar, const Vector &etax, const Vector &xix, Vector *result) const
    {
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::ScalarVectorAddVector has not been overridden!\n");
        #endif
        *result = xix;
        result->AlphaXaddThis(scalar, etax);
        return *result;
    };

    Vector &Manifold::VectorLinearCombination(const Variable &x, realdp scalar1, const Vector &etax, realdp scalar2, const Vector &xix, Vector *result) const
    {
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::VectorLinearCombination has not been overridden!\n");
        #endif
        *result = xix;
        result->ScalarTimesThis(scalar2);
        result->AlphaXaddThis(scalar1, etax);
        
        return *result;
    };

    Variable &Manifold::Retraction(const Variable &x, const Vector &etax, Variable *result) const
    {
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::Retraction has not been overridden!\n");
        #endif
        *result = x; result->AlphaXaddThis(1, etax); /*result = x + etax*/
        return *result;
    };

    Vector &Manifold::InvRetraction(const Variable &x, const Variable &y, Vector *result) const
    {
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::InvRetraction has not been overridden!\n");
        #endif
        *result = y; result->AlphaXaddThis(-1, x); /*result = y - x*/
        
        return ExtrProjection(x, *result, result);
    };

	Vector &Manifold::coTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
	{
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::coTangentVector has not been overridden!\n");
        #endif
        
        return Projection(x, xiy, result);
	};

	Vector &Manifold::DiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir) const
	{
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::DiffRetraction has not been overridden!\n");
        #endif
        
        Projection(y, xix, result);

        if (IsEtaXiSameDir && HasHHR)
        {
            Vector beta(3);
            realdp *betaptr = beta.ObtainWriteEntireData();
            realdp EtatoXi = sqrt(Metric(x, etax, etax) / Metric(x, xix, xix));
            betaptr[0] = std::sqrt(Metric(x, etax, etax) / Metric(x, *result, *result)) / EtatoXi;
            betaptr[1] = Metric(x, etax, etax);
            betaptr[2] = Metric(x, *result, *result) * EtatoXi * EtatoXi;
            etax.AddToFields("beta", beta);
            
            if (HasHHR)
            {
                etax.AddToFields("betaTReta", (*result) * (betaptr[0] * EtatoXi));
            }
        }
        
        return *result;
	};

	realdp Manifold::Beta(const Variable &x, const Vector &etax) const
	{
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::Beta has not been overridden!\n");
        #endif
		return 1;
	};

	realdp Manifold::Dist(const Variable &x1, const Variable &x2) const
	{
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::Dist has not been overridden!\n");
        #endif
		return (x1 - x2).Fnorm();
	};

	Vector &Manifold::VectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const
	{
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::VectorTransport has not been overridden!\n");
        #endif
        if (HasHHR)
            return LCVectorTransport(x, etax, y, xix, result);
        
        return Projection(y, xix, result);
	};

	Vector &Manifold::InverseVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
	{
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::InverseVectorTransport has not been overridden!\n");
        #endif
		if (HasHHR)
			return LCInverseVectorTransport(x, etax, y, xiy, result);
		
		return Projection(x, xiy, result);
	};

	LinearOPE &Manifold::HInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, integer start, integer end, LinearOPE *result) const
	{
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::HInvTran has not been overridden!\n");
        #endif
		if (HasHHR)
			return LCHInvTran(x, etax, y, Hx, start, end, result);
		
        *result = Hx;
        return *result;
	};

	LinearOPE &Manifold::TranH(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, integer start, integer end, LinearOPE *result) const
	{
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::TranH has not been overridden!\n");
        #endif
		if (HasHHR)
			return LCTranH(x, etax, y, Hx, start, end, result);
		
        *result = Hx;
		return *result;
	};

	LinearOPE &Manifold::TranHInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, LinearOPE *result) const
	{
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::TranHInvTran has not been overridden!\n");
        #endif
        
		if (HasHHR)
			return LCTranHInvTran(x, etax, y, Hx, result);
        
        *result = Hx;
        return *result;
	};

	LinearOPE &Manifold::HaddScaledRank1OPE(const Variable &x, const LinearOPE &Hx, realdp scalar, const Vector &etax, const Vector &xix, LinearOPE *result) const
	{
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::HaddScaledRank1OPE has not been overridden!\n");
        #endif
        /*even if etax and xix are complex elements, the real reparameterizations are used. In other words, we viewed it is
        2 n dimension real data. Therefore, the low rank update uses real operations. */
        bool etaxiscomplex = etax.Getiscomplex(), xixiscomplex = xix.Getiscomplex();
        etax.Setiscomplex(false);
        xix.Setiscomplex(false);
        *result = Hx;
        Vector etaxreshape = etax; etaxreshape.Reshape(etax.Getlength());
        Vector xixreshape = xix; xixreshape.Reshape(xix.Getlength());
        result->HaddRankone(scalar, etaxreshape, xixreshape);
        etax.Setiscomplex(etaxiscomplex);
        xix.Setiscomplex(xixiscomplex);
        return *result;
	};

	Vector &Manifold::ObtainEtaxFlat(const Variable &x, const Vector &etax, Vector *result) const
	{
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::ObtainEtaxFlat has not been overridden!\n");
        #endif
        *result = etax;
        return *result;
	};

	Vector &Manifold::ObtainIntr(const Variable &x, const Vector &etax, Vector *result) const
	{
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::ObtainIntr has not been overridden!\n");
        #endif
        *result = etax;
        return *result;
	};

	Vector &Manifold::ObtainExtr(const Variable &x, const Vector &intretax, Vector *result) const
	{
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::ObtainExtr has not been overridden!\n");
        #endif
        *result = intretax;
        return *result;
	};

	Vector &Manifold::ObtainNorVerIntr(const Variable &x, const Vector &etax, Vector *result) const
	{
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::ObtainNorVerIntr has not been overridden!\n");
        #endif
        *result = etax;
        return *result;
	};

	Vector &Manifold::ObtainNorVerExtr(const Variable &x, const Vector &intretax, Vector *result) const
	{
        #ifdef CHECKMANIFOLDOVERRIDDEN
                printf("Manifold::ObtainNorVerExtr has not been overridden!\n");
        #endif
        *result = intretax;
        return *result;
	};

    void Manifold::CheckParams(void) const
    {
        printf("GENERAL PARAMETERS:\n");
        printf("name          :%15s,\n", name.c_str());
        printf("IsIntrApproach:%15d,\t", IsIntrApproach);
        printf("IntrinsicDim  :%15d,\n", IntrinsicDim);
        printf("ExtrinsicDim  :%15d,\t", ExtrinsicDim);
        printf("HasHHR        :%15d,\n", HasHHR);
    };

    void Manifold::SetParams(PARAMSMAP params)
    {
        PARAMSMAP::iterator iter;
        for (iter = params.begin(); iter != params.end(); iter++)
        {
            if (iter->first == static_cast<std::string> ("HasHHR"))
            {
                SetHasHHR(((static_cast<integer> (iter->second)) != 0));
            }
        }
    };

	void Manifold::CheckIntrExtr(Variable x) const
	{
		printf("==============Check Intrinsic/Extrinsic transform=========\n");
		Vector exetax(EMPTYEXTR);
		Vector inetax(EMPTYINTR);

		x.Print("x");
		exetax.RandGaussian();
		ExtrProjection(x, exetax, &exetax);
		exetax.Print("exetax1");
		 ObtainIntr(x, exetax, &inetax);
        bool Isintr = GetIsIntrinsic();
        SetIsIntrApproach(false);
		printf("extr inp:%g\n", Metric(x, exetax, exetax));
        SetIsIntrApproach(true);
		printf("intr inp:%g\n", Metric(x, inetax, inetax));
        SetIsIntrApproach(Isintr);
		inetax.Print("inetax1");
		ObtainExtr(x, inetax, &exetax);
		exetax.Print("exetax2");
		ObtainIntr(x, exetax, &inetax);
		inetax.Print("inetax2");
		printf("exeta1 and inetax1 should approximately equal exetax2 and inetax2 respectively!\n");
	};

	void Manifold::CheckRetraction(Variable x) const
	{
		printf("==============Check Retraction=========\n");
		Vector etax(EMPTYEXTR), FDetax(EMPTYEXTR);
		etax.RandGaussian();
		ExtrProjection(x, etax, &etax);
        ScalarTimesVector(x, 1.0 / sqrt(Metric(x, etax, etax)), etax, &etax);
		x.Print("x:");
		etax.Print("etax:");
#ifdef SINGLE_PRECISION
		realdp eps = static_cast<realdp> (1e-3);
#else
        realdp eps = static_cast<realdp> (1e-5);
#endif
		Variable y(x);
		ScalarTimesVector(x, eps, etax, &etax);
		if (IsIntrApproach)
		{
			Vector inetax(EMPTYINTR);
			ObtainIntr(x, etax, &inetax);
			Retraction(x, inetax, &y);
		}
		else
		{
			Retraction(x, etax, &y);
		}
        FDetax = (y - x) / eps;
		FDetax.Print("FDetax:");

		printf("etax should approximately equal FDetax = (R(eps etax)-R(etax))/eps!\n");
	};

	void Manifold::CheckDiffRetraction(Variable x, bool IsEtaXiSameDir) const
	{
		printf("==============Check Differentiated Retraction=========\n");
		Vector etax(EMPTYEXTR), xix(EMPTYEXTR), zetax(EMPTYEXTR);
		etax.RandGaussian();
		ExtrProjection(x, etax, &etax);
		if (IsEtaXiSameDir)
		{
            xix = etax;
		}
		else
		{
			xix.RandGaussian();
			ExtrProjection(x, xix, &xix);
		}
		x.Print("x:");
		etax.Print("etax:");
		Variable y(x);
		if (IsIntrApproach)
		{
			Vector inetax(EMPTYINTR), inxix(EMPTYINTR), inzetax(EMPTYINTR);
			ObtainIntr(x, etax, &inetax);
			ObtainIntr(x, xix, &inxix);
            std::cout << "************************************************************************************************************" << std::endl;
			Retraction(x, inetax, &y);
			DiffRetraction(x, inetax, y, inxix, &inzetax, IsEtaXiSameDir);
			ObtainExtr(y, inzetax, &zetax);
		}
		else
		{
			Retraction(x, etax, &y);
			DiffRetraction(x, etax, y, xix, &zetax, IsEtaXiSameDir);
		}
		y.Print("y:");
		zetax.Print("zetax:");
		Variable yeps(x);
		realdp eps = static_cast<realdp> (1e-5);
        ScalarVectorAddVector(x, eps, xix, etax, &etax);
		if (IsIntrApproach)
		{
			Vector inetax(EMPTYINTR);
			ObtainIntr(x, etax, &inetax);
			Retraction(x, inetax, &yeps);
		}
		else
		{
			Retraction(x, etax, &yeps);
		}
        zetax = (yeps - y) / eps;
		ExtrProjection(y, zetax, &zetax);
		zetax.Print("FDzetax:");
		printf("zetax = T_{R_etax} xix should approximately equal FDzetax = (R(etax+eps xix) - R(etax))/eps!\n");
	};

	void Manifold::CheckLockingCondition(Variable x) const
	{
		printf("==============Check Locking Condition=========\n");
		Vector etax(EMPTYEXTR), xix(EMPTYEXTR), zetax(EMPTYEXTR);
		etax.RandGaussian();
		ExtrProjection(x, etax, &etax);
		ScalarTimesVector(x, 1, etax, &xix); /* genrandreal() + static_cast<realdp> (0.5) */
		Variable y(x);
		if (IsIntrApproach)
		{
			Vector inetax(EMPTYINTR), inxix(EMPTYINTR), inzetax(EMPTYINTR);
			ObtainIntr(x, etax, &inetax);
			ObtainIntr(x, xix, &inxix);
			Retraction(x, inetax, &y);
			DiffRetraction(x, inetax, y, inxix, &inzetax, true);
			if (inetax.FieldsExist("beta"))
			{
                Vector beta = inetax.Field("beta");
                const realdp *betav = beta.ObtainReadData();
                printf("beta = |etax| / |T_{etax} etax|: %g\n", betav[0]);
			}
			else
			{
				printf("beta: %d\n", 1);
			}
			printf("|xix| / |T_{etax} xix|:%g\n", sqrt(Metric(x, inxix, inxix) / Metric(x, inzetax, inzetax)));
			ScalarTimesVector(x, sqrt(Metric(x, inxix, inxix) / Metric(x, inzetax, inzetax)), inzetax, &inzetax);
			ObtainExtr(y, inzetax, &zetax);
			zetax.Print("Beta DiffRetraction zetax:");
			VectorTransport(x, inetax, y, inxix, &inzetax);
            
			ObtainExtr(y, inzetax, &zetax);
			zetax.Print("Vector Transport zetax:");
		}
		else
		{
			Retraction(x, etax, &y);
			DiffRetraction(x, etax, y, xix, &zetax, true);
			if (etax.FieldsExist("beta"))
			{
                Element beta = etax.Field("beta");
                const realdp *betav = beta.ObtainReadData();
                printf("beta = |etax| / |T_{etax} etax|: %g\n", betav[0]);
			}
			else
			{
				printf("beta: %d\n", 1);
			}
			printf("|xix| / |T_{etax} xix|:%g\n", sqrt(Metric(x, xix, xix) / Metric(y, zetax, zetax)));
			ScalarTimesVector(y, sqrt(Metric(x, xix, xix) / Metric(y, zetax, zetax)), zetax, &zetax);
			zetax.Print("Beta DiffRetraction zetax:");
			VectorTransport(x, etax, y, xix, &zetax);
			zetax.Print("Vector Transport zetax:");
		}
		printf("Beta DiffRetraction zetax should approximately equal Vector Transport zetax!\n");
	};

	void Manifold::CheckcoTangentVector(Variable x) const
	{
		printf("==============Check CoTangentVector=========\n");
		Vector etax(EMPTYEXTR), xix(EMPTYEXTR), zetay(EMPTYEXTR), xiy(EMPTYEXTR), zetax(EMPTYEXTR);
		etax.RandGaussian();
		ExtrProjection(x, etax, &etax);

		xix.RandGaussian();
		ExtrProjection(x, xix, &xix);

		Variable y(x);
		if (IsIntrApproach)
		{
			Vector inetax(EMPTYINTR), inxix(EMPTYINTR), inzetay(EMPTYINTR), inxiy(EMPTYINTR), inzetax(EMPTYINTR);
            ObtainIntr(x, etax, &inetax);
			ObtainIntr(x, xix, &inxix);
			Retraction(x, inetax, &y);
			DiffRetraction(x, inetax, y, inxix, &inzetay, false);
			ObtainExtr(y, inzetay, &zetay);

			xiy.RandGaussian();
            ExtrProjection(y, xiy, &xiy);
			ObtainIntr(y, xiy, &inxiy);
			printf("<xiy, T_{R_{eta}} xix>:%g\n", Metric(y, inxiy, inzetay));

			coTangentVector(x, inetax, y, inxiy, &inzetax);
			ObtainExtr(x, inzetax, &zetax);
			printf("C(x, etax, xiy) [xix]:%g\n", Metric(x, inzetax, inxix));
		}
		else
		{
			Retraction(x, etax, &y);
			DiffRetraction(x, etax, y, xix, &zetay, false);
			xiy.RandGaussian();
			ExtrProjection(y, xiy, &xiy);
			ScalarTimesVector(y, sqrt(Metric(y, xiy, xiy)), xiy, &xiy);
			printf("<xiy, T_{R_{eta}} xix>:%g\n", Metric(y, xiy, zetay));
			coTangentVector(x, etax, y, xiy, &zetax);
			printf("C(x, etax, xiy) [xix]:%g\n", Metric(x, zetax, xix));
		}
		printf("<xiy, T_{R_{eta}} xix> should approximately equal C(x, etax, xiy) [xix]!\n");
	};

	void Manifold::CheckIsometryofVectorTransport(Variable x) const
	{
		printf("==============Check Isometry of the Vector Transport=========\n");
		Vector etax(EMPTYEXTR), xix(EMPTYEXTR), zetay(EMPTYEXTR);
		etax.RandGaussian();
		ExtrProjection(x, etax, &etax);

		xix.RandGaussian();
		ExtrProjection(x, xix, &xix);

		Variable y(x);
		if (IsIntrApproach)
		{
			Vector inetax(EMPTYINTR), inxix(EMPTYINTR), inzetay(EMPTYINTR);
			ObtainIntr(x, etax, &inetax);
			ObtainIntr(x, xix, &inxix);
			Retraction(x, inetax, &y);
			VectorTransport(x, inetax, y, inxix, &inzetay);
			printf("Before vector transport:%g, After vector transport:%g\n", Metric(x, inxix, inxix), Metric(y, inzetay, inzetay));
		}
		else
		{
			Retraction(x, etax, &y);
			VectorTransport(x, etax, y, xix, &zetay);
			y.Print("y:");
			zetay.Print("zetay:");
			printf("Before vector transport:%g, After vector transport:%g\n", Metric(x, xix, xix), Metric(y, zetay, zetay));
		}
		printf("|xix| (Before vector transport) should approximately equal |T_{R_etax} xix| (After vector transport)\n");
	};

	void Manifold::CheckIsometryofInvVectorTransport(Variable x) const
	{
		printf("==============Check Isometry of the Inverse Vector Transport=========\n");
		Vector etax(EMPTYEXTR), xix(EMPTYEXTR), zetay(EMPTYEXTR);

		etax.RandGaussian();
		ExtrProjection(x, etax, &etax);

		Variable y(x);
		if (IsIntrApproach)
		{
			Vector inetax(EMPTYINTR), inxix(EMPTYINTR), inzetay(EMPTYINTR);
			ObtainIntr(x, etax, &inetax);
			Retraction(x, inetax, &y);
			zetay.RandGaussian();
			ExtrProjection(y, zetay, &zetay);
			ScalarTimesVector(y, sqrt(Metric(y, zetay, zetay)), zetay, &zetay);
			ObtainIntr(y, zetay, &inzetay);

			InverseVectorTransport(x, inetax, y, inzetay, &inxix);
			printf("Before inverse vector transport:%g, After inverse vector transport:%g\n", Metric(y, inzetay, inzetay), Metric(x, inxix, inxix));
		}
		else
		{
			Retraction(x, etax, &y);
			zetay.RandGaussian();
			ExtrProjection(y, zetay, &zetay);
			InverseVectorTransport(x, etax, y, zetay, &xix);
			x.Print("x:");
			xix.Print("xix:");
			printf("Before inverse vector transport:%g, After inverse vector transport:%g\n", Metric(y, zetay, zetay), Metric(x, xix, xix));
		}
		printf("|zetay| (Before inverse vector transport) should approximately equal |T_{R_etax}^{-1} zetay| (After inverse vector transport)\n");
	};

	void Manifold::CheckVecTranComposeInverseVecTran(Variable x) const
	{
		printf("==============Check Vector Transport Compose Inverse Vector Transport=========\n");
		Vector etax(EMPTYEXTR), xix(EMPTYEXTR), zetay(EMPTYEXTR);

		etax.RandGaussian();
		ExtrProjection(x, etax, &etax);
		xix.RandGaussian();
		ExtrProjection(x, xix, &xix);

		Variable y(x);
		if (IsIntrApproach)
		{
			Vector inetax(EMPTYINTR), inxix(EMPTYINTR), inzetay(EMPTYINTR);
			ObtainIntr(x, etax, &inetax);
			Retraction(x, inetax, &y);
			ObtainIntr(x, xix, &inxix);
			xix.Print("xix:");
			VectorTransport(x, inetax, y, inxix, &inzetay);
			InverseVectorTransport(x, inetax, y, inzetay, &inxix);
			ObtainExtr(x, inxix, &xix);
			xix.Print("T^{-1} ciric T xix:");
			printf("xix and T^{-1} ciric T xix should be similar!\n");
		}
		else
		{
			Retraction(x, etax, &y);
			xix.Print("xix:");
			VectorTransport(x, etax, y, xix, &zetay);
			InverseVectorTransport(x, etax, y, zetay, &xix);
			xix.Print("T^{-1} ciric T xix:");
			printf("xix and T^{-1} ciric T xix should be similar!\n");
		}
	};

	void Manifold::CheckTranHInvTran(Variable x) const
	{
		printf("==============Check Transport of a Hessian approximation=========\n");
		Vector etax(EMPTYEXTR);
		Variable y(x);

		etax.RandGaussian();
		ExtrProjection(x, etax, &etax);

		if (IsIntrApproach)
		{
            LinearOPE Hx(EMPTYINTR.Getlength(), EMPTYINTR.Getlength()), result(EMPTYINTR.Getlength(), EMPTYINTR.Getlength());
			Vector inetax(EMPTYINTR);
			ObtainIntr(x, etax, &inetax);
			Retraction(x, inetax, &y);
			Hx.ScaledIdOPE();
			Hx.Print("Hx before:");
			TranHInvTran(x, inetax, y, Hx, &result);
			result.Print("Hx after:");
		}
		else
		{
            LinearOPE Hx(EMPTYEXTR.Getlength(), EMPTYEXTR.Getlength()), result(EMPTYEXTR.Getlength(), EMPTYEXTR.Getlength());
			Hx.ScaledIdOPE();
			Hx.Print("Hx before:");
			Retraction(x, etax, &y);
			Vector zetay1(EMPTYEXTR), zetay2(EMPTYEXTR);
			zetay1.RandGaussian();
			ExtrProjection(y, zetay1, &zetay1);
			TranHInvTran(x, etax, y, Hx, &result);
			result.Print("Hx after:");
			zetay1.Print("zetay:");
            LinearOPEEta(y, result, zetay1, &zetay2);
			zetay2.Print("Hx zetay:");
		}
	};

	void Manifold::CheckHaddScaledRank1OPE(Variable x) const
	{
		printf("==============Check Rank one Update to a Hessian Approximation=========\n");
		realdp scalar = 1.0;
		Vector etax(EMPTYEXTR), xix(EMPTYEXTR);
		etax.RandGaussian();
		ExtrProjection(x, etax, &etax);
		xix.RandGaussian();
		ExtrProjection(x, xix, &xix);

		if (IsIntrApproach)
		{
            LinearOPE Hx(EMPTYINTR.Getlength(), EMPTYINTR.Getlength()), result(EMPTYINTR.Getlength(), EMPTYINTR.Getlength());
			Vector inetax(EMPTYINTR), inxix(EMPTYINTR);
			ObtainIntr(x, etax, &inetax);
			ObtainIntr(x, xix, &inxix);
			Hx.ScaledIdOPE();
			Hx.Print("Hx before:");
			HaddScaledRank1OPE(x, Hx, scalar, inetax, inxix, &result);
			inetax.Print("etax:");
			inxix.Print("xix:");
			result.Print("Hx after:");
		}
		else
		{
            LinearOPE Hx(EMPTYEXTR.Getlength(), EMPTYEXTR.Getlength()), result(EMPTYEXTR.Getlength(), EMPTYEXTR.Getlength());
			Hx.ScaledIdOPE();
			Hx.Print("Hx before:");
			HaddScaledRank1OPE(x, Hx, scalar, etax, xix, &result);
			etax.Print("etax:");
			xix.Print("xix:");
			result.Print("Hx after:");
		}
	};

    void Manifold::Obtainnu1nu2forLC(const Variable &x, const Vector &etax, const Variable &y) const
    {
        Vector eps1(etax), nu1(etax), nu2(etax);
        
        if (!etax.FieldsExist("beta") || !etax.FieldsExist("betaTReta"))
        {
            DiffRetraction(x, etax, y, etax, &eps1, true);
        }
        Vector TRetaVector = etax.Field("betaTReta");
        HasHHR = false; VectorTransport(x, etax, y, etax, &eps1); HasHHR = true;
        ScalarTimesVector(y, sqrt(Metric(x, etax, etax) / Metric(y, eps1, eps1)), eps1, &eps1); /*note \|etax\|_x = \|TRetaVector\|_y*/
        Element tau1tau2(2);
        realdp *tau1tau2ptr = tau1tau2.ObtainWriteEntireData();
        ScalarTimesVector(y, 2.0, eps1, &nu1);
        VectorLinearCombination(y, -1.0, eps1, -1.0, TRetaVector, &nu2);
        tau1tau2ptr[0] = static_cast<realdp> (2) / Metric(y, nu1, nu1);
        tau1tau2ptr[1] = static_cast<realdp> (2) / Metric(y, nu2, nu2);
        
        etax.AddToFields("tau1tau2", tau1tau2);
        etax.AddToFields("nu1", nu1);
        etax.AddToFields("nu2", nu2);
    };

    Vector &Manifold::LCVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const
    {
        if (!etax.FieldsExist("tau1tau2"))
        {
            Obtainnu1nu2forLC(x, etax, y);
        }
        HasHHR = false; VectorTransport(x, etax, y, xix, result); HasHHR = true;
        ScalarTimesVector(y, sqrt(Metric(x, xix, xix) / Metric(y, *result, *result)), *result, result);
        
        Element tau1tau2 = etax.Field("tau1tau2");
        const realdp *tau1tau2ptr = tau1tau2.ObtainReadData();
        Vector nu1 = etax.Field("nu1");
        Vector nu2 = etax.Field("nu2");
        realdp temp = - Metric(y, *result, nu1);
        ScalarVectorAddVector(y, temp * tau1tau2ptr[0], nu1, *result, result);
        temp = -Metric(y, *result, nu2);
        ScalarVectorAddVector(y, temp * tau1tau2ptr[1], nu2, *result, result);
        
        return *result;
    };

    Vector &Manifold::LCInverseVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
    {
        if (!etax.FieldsExist("tau1tau2"))
        {
            Obtainnu1nu2forLC(x, etax, y);
        }
        
        Element tau1tau2 = etax.Field("tau1tau2");
        const realdp *tau1tau2ptr = tau1tau2.ObtainReadData();
        Vector nu1 = etax.Field("nu1");
        Vector nu2 = etax.Field("nu2");
        realdp temp = -Metric(y, xiy, nu2);
        VectorLinearCombination(y, temp * tau1tau2ptr[1], nu2, 1, xiy, result);
        temp = -Metric(y, *result, nu1);
        ScalarVectorAddVector(y, temp * tau1tau2ptr[0], nu1, *result, result);
        
        HasHHR = false; InverseVectorTransport(x, etax, y, *result, result); HasHHR = true;
        return *result;
    };

    LinearOPE &Manifold::LCHInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, integer start, integer end, LinearOPE *result) const
    {
        if (!etax.FieldsExist("tau1tau2"))
        {
            Obtainnu1nu2forLC(x, etax, y);
        }
        Element tau1tau2 = etax.Field("tau1tau2");
        const realdp *tau1tau2ptr = tau1tau2.ObtainReadData();
        Vector nu1 = etax.Field("nu1");
        Vector nu2 = etax.Field("nu2");
        const realdp *nu1TV = nu1.ObtainReadData();
        const realdp *nu2TV = nu2.ObtainReadData();
        
        HasHHR = false; HInvTran(x, etax, y, Hx, start, end, result); HasHHR = true;
        
        realdp *resultTV = result->ObtainWritePartialData();
        char *sider = const_cast<char *> ("r");
        integer ell = Hx.Getsize()[0], length = etax.Getlength();
        realdp *work = new realdp[ell];

        /* resultTV(:, start : start + length - 1) <- resultTV(:, start : start + length - 1) * (I - tau1tau2(0) * nu1TV * nu1TV^T),
        details: www.netlib.org/lapack/explore-html/db/d10/larfx_8f.html */
        larfx_(sider, &ell, &length, const_cast<realdp *> (nu1TV), const_cast<realdp *> (tau1tau2ptr), resultTV + start * ell, &ell, work);
        /* resultTV(:, start : start + length - 1) <- resultTV(:, start : start + length - 1) * (I - tau1tau2(1) * nu2TV * nu2TV^T),
        details: www.netlib.org/lapack/explore-html/db/d10/larfx_8f.html */
        larfx_(sider, &ell, &length, const_cast<realdp *> (nu2TV), const_cast<realdp *> (tau1tau2ptr + 1), resultTV + start * ell, &ell, work);
        delete[] work;
        return *result;
    };

    LinearOPE &Manifold::LCTranH(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, integer start, integer end, LinearOPE *result) const
    {
        if (!etax.FieldsExist("tau1tau2"))
        {
            Obtainnu1nu2forLC(x, etax, y);
        }
        Element tau1tau2 = etax.Field("tau1tau2");
        const realdp *tau1tau2ptr = tau1tau2.ObtainReadData();
        Vector nu1 = etax.Field("nu1");
        Vector nu2 = etax.Field("nu2");
        const realdp *nu1TV = nu1.ObtainReadData();
        const realdp *nu2TV = nu2.ObtainReadData();
        HasHHR = false; TranH(x, etax, y, Hx, start, end, result); HasHHR = true;
        realdp *resultTV = result->ObtainWritePartialData();
        
        char *sidel = const_cast<char *> ("l");
        integer ell = Hx.Getsize()[0], length = etax.Getlength();
        realdp *work = new realdp[ell];
        /* resultTV(start : start + length - 1, :) <- (I - tau1tau2(0) * nu1TV * nu1TV^T) * resultTV(start : start + length - 1, :),
        details: www.netlib.org/lapack/explore-html/db/d10/larfx_8f.html */
        larfx_(sidel, &length, &ell, const_cast<realdp *> (nu1TV), const_cast<realdp *> (tau1tau2ptr), resultTV + start, &ell, work);
        /* resultTV(start : start + length - 1, :) <- (I - tau1tau2(1) * nu2TV * nu2TV^T) * resultTV(start : start + length - 1, :),
        details: www.netlib.org/lapack/explore-html/db/d10/larfx_8f.html */
        larfx_(sidel, &length, &ell, const_cast<realdp *> (nu2TV), const_cast<realdp *> (tau1tau2ptr + 1), resultTV + start, &ell, work);
        delete[] work;
        return *result;
    };

    LinearOPE &Manifold::LCTranHInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, LinearOPE *result) const
    {
        if (!etax.FieldsExist("tau1tau2"))
        {
            Obtainnu1nu2forLC(x, etax, y);
        }
        
        Element tau1tau2 = etax.Field("tau1tau2");
        const realdp *tau1tau2ptr = tau1tau2.ObtainReadData();
        Vector nu1 = etax.Field("nu1");
        Vector nu2 = etax.Field("nu2");
        const realdp *nu1TV = nu1.ObtainReadData();
        const realdp *nu2TV = nu2.ObtainReadData();
        HasHHR = false; TranHInvTran(x, etax, y, Hx, result); HasHHR = true;
        realdp *resultTV = result->ObtainWritePartialData();
        
        char *sidel = const_cast<char *> ("l"), *sider = const_cast<char *> ("r");
        integer ell = Hx.Getsize()[0], length = etax.Getlength();
        realdp *work = new realdp[ell];
        /* resultTV <- resultTV * (I - tau1tau2(0) * nu1TV * nu1TV^T),
        details: www.netlib.org/lapack/explore-html/db/d10/larfx_8f.html */
        larfx_(sider, &ell, &length, const_cast<realdp *> (nu1TV), const_cast<realdp *> (tau1tau2ptr), resultTV, &ell, work);
        /* resultTV <- resultTV * (I - tau1tau2(1) * nu2TV * nu2TV^T),
        details: www.netlib.org/lapack/explore-html/db/d10/larfx_8f.html */
        larfx_(sider, &ell, &length, const_cast<realdp *> (nu2TV), const_cast<realdp *> (tau1tau2ptr + 1), resultTV, &ell, work);
        /* resultTV <- (I - tau1tau2(0) * nu1TV * nu1TV^T) * resultTV,
        details: www.netlib.org/lapack/explore-html/db/d10/larfx_8f.html */
        larfx_(sidel, &length, &ell, const_cast<realdp *> (nu1TV), const_cast<realdp *> (tau1tau2ptr), resultTV, &ell, work);
        /* resultTV <- (I - tau1tau2(1) * nu2TV * nu2TV^T) * resultTV,
        details: www.netlib.org/lapack/explore-html/db/d10/larfx_8f.html */
        larfx_(sidel, &length, &ell, const_cast<realdp *> (nu2TV), const_cast<realdp *> (tau1tau2ptr + 1), resultTV, &ell, work);
        delete[] work;
        return *result;
    };

    Vector &Manifold::TangentSpaceProximalMap(Variable &x, const Vector &etax, realdp adavalue, realdp SMtol, realdp SMlambda, const Problem *prob, Vector *inoutinitD, integer *outSMiter, integer *outSMCGiter, Vector *result) const
    { /*eta = argmin g(gfx, eta) + 0.5 \|eta\|_W^2 + h(eta), W is either a scalar or a weight matrix that is given by function "PreConditioner" defined in Problem.
              Only extrinsic representation is supported. */
        integer dimNorVec = ExtrinsicDim - IntrinsicDim;
        if(inoutinitD->GetSpace() == nullptr)
        {
            *inoutinitD = Vector (dimNorVec);
            inoutinitD->SetToZeros();
        }
        Vector SMz(dimNorVec), Fz(dimNorVec), SMd(dimNorVec), SMu(dimNorVec), Fu(dimNorVec);
        
        /*parameters*/
        realdp SMtau = 0.1, SMalpha = 0.1;
        integer SMmaxiter = 100, maxSMbtiter = 10;
        std::string status;
        integer SMbtiter = 0;
        bool earlytermination = false;
        Vector Weight;
        prob->PreConditioner(x, Weight, &Weight);
        Weight.ScalarTimesThis(adavalue); /*Weight = adavalue * Weight*/
        
        realdp nFz = 0, nFu = 0, mu = 0;
        SMz = *inoutinitD;
        Vector BLambda(x);
        ComputeBLambda(x, SMz, Weight, etax, prob, &BLambda);
        
        EW(x, BLambda, Weight, prob, &Fz);
        nFz = Fz.Fnorm();
        
        SMd.SetToZeros();
        
        (*outSMiter) = 0;
        (*outSMCGiter) = 0;
        integer outSMCGiter_i = 0;
        while(nFz * nFz > SMtol && (*outSMiter) < SMmaxiter)
        {
            mu = ((nFz < 0.1) ? nFz : 0.1);
            mu = ((mu > 1e-11) ? mu : 1e-11);
            mu = SMlambda *  mu;
            myCG(x, Fz, dimNorVec, mu, SMtau, SMlambda * nFz, ((dimNorVec < 30) ? dimNorVec : 30), Weight, BLambda, SMd, prob, &outSMCGiter_i, &SMd);
            (*outSMCGiter) = (*outSMCGiter) + outSMCGiter_i;
            SMu = SMz; SMu.AlphaXaddThis(1, SMd); /* SMu = SMz + SMd; */
            ComputeBLambda(x, SMu, Weight, etax, prob, &BLambda);
            EW(x, BLambda, Weight, prob, &Fu);
            nFu = Fu.Fnorm();

            SMalpha = 1;
            SMbtiter = 0;
            while(nFu * nFu > nFz * nFz * (1 - 0.001 * SMalpha) && SMbtiter < maxSMbtiter)
            {
                SMalpha *= 0.5;
                SMu = SMz; SMu.AlphaXaddThis(SMalpha, SMd); /* SMu = SMalpha * SMd + SMz; */
                ComputeBLambda(x, SMu, Weight, etax, prob, &BLambda);
                EW(x, BLambda, Weight, prob, &Fu);
                nFu = Fu.Fnorm();
                SMbtiter++;
            }
            /*We let the SSN run for at least 1 iteration before early termination*/
            if(SMbtiter == maxSMbtiter && (*outSMiter) != 0)
            {
                earlytermination = true;
                break;
            }
            SMz = SMu;
            Fz = Fu;
            nFz = nFu;
            (*outSMiter)++;
        }
        *inoutinitD = SMz;
        prob->ProxW(BLambda, Weight, result); result->AlphaXaddThis(-1, x); /*result = prob->ProxW(BLambda, Weight) - x*/
        Projection(x, *result, result);
        return *result;
    };

    Vector &Manifold::myCG(const Variable &x, const Vector &nb, integer dimNorVec, realdp mu, realdp tau, realdp lambdanFz, integer maxiter, const Vector &Weight, const Vector &BLambda, const Vector &init, const Problem *prob, integer *CGiter, Vector *result) const
    {
        Vector r(dimNorVec), p(dimNorVec), Ap(dimNorVec), SMd(init);
        realdp rr0, alpha, beta;
        
        GLdW(x, SMd, Weight, BLambda, prob, &Ap).AlphaXaddThis(mu, SMd); /* Ap = mu * SMd + GLdW(x, SMd, Weight, BLambda, prob, &tmp); */
        *result = SMd;
        r = Ap; r.ScalarTimesThis(-1); r.AlphaXaddThis(-1, nb); /*r = -1 * nb - Ap; r is negative residual*/
        p = r;
        (*CGiter) = 0;
        while(r.Fnorm() > tau * ((lambdanFz * result->Fnorm() < 1.0) ? lambdanFz * result->Fnorm() : 1.0) && (*CGiter) < maxiter)
        {
            GLdW(x, p, Weight, BLambda, prob, &Ap).AlphaXaddThis(mu, p); /* Ap = mu * p + GLdW(x, p, Weight, BLambda, prob, &tmp); */
            rr0 = r.DotProduct(r);
            alpha = rr0 / p.DotProduct(Ap);
            result->AlphaXaddThis(alpha, p); /* (*result) = alpha * p + (*result); */
            r.AlphaXaddThis(- alpha, Ap); /* r = r - alpha * Ap; */
            beta = r.DotProduct(r) / rr0;
            p.ScalarTimesThis(beta); p.AlphaXaddThis(1, r); /* p = r + beta * p; */
            (*CGiter)++;
        }
        return *result;
    };

    Vector &Manifold::ComputeBLambda(const Variable &x, const Vector &d, const Vector &Weight, const Vector &gfx, const Problem *prob, Vector *result) const
    {
        calAadj(x, d, prob, result);
        result->ScalarTimesThis(-1); result->AlphaXaddThis(1, gfx); /*result = gfx - calAadj(x, d, prob) */
        
        integer nblock = Weight.Getlength();
        integer blocksize = x.Getlength() / nblock;
        realdp *resultptr = result->ObtainWritePartialData();
        
        for(integer i = 0; i < nblock; i++)
        {
            const realdp L = Weight.ObtainReadData()[i];
            realdp tmp = 1.0 / L;
            scal_(&blocksize, &tmp, resultptr + i * blocksize, &GLOBAL::IONE);
        }
        result->ScalarTimesThis(-1); result->AlphaXaddThis(1, x); /* result = x - result; */
        return *result;
    };

    Vector &Manifold::EW(const Variable &x, const Vector &BLambda, const Vector &Weight, const Problem *prob, Vector *result) const
    {
        Vector DLambda(BLambda);
        prob->ProxW(BLambda, Weight, &DLambda);
        DLambda.AlphaXaddThis(-1, x); /* DLambda = prob->ProxW(BLambda, Weight) - x; */
        return calA(x, DLambda, prob, result);
    };

    Vector &Manifold::GLdW(const Variable &x, const Vector &d, const Vector &Weight, const Vector &BLambda, const Problem *prob, Vector *result) const
    {
        Vector VecTMP1(x), VecTMP2(x);
        calAadj(x, d, prob, &VecTMP1);
        prob->CalJW(BLambda, VecTMP1, Weight, &VecTMP2);
        
        integer nblock = Weight.Getlength();
        integer blocksize = x.Getlength() / nblock;
        realdp *Vec2ptr = VecTMP2.ObtainWritePartialData();
        
        for(integer i = 0; i < nblock; i++)
        {
            const realdp L = Weight.ObtainReadData()[i];
            realdp tmp = 1.0 / L;
            scal_(&blocksize, &tmp, Vec2ptr + i * blocksize, &GLOBAL::IONE);
        }
        
        return calA(x, VecTMP2, prob, result);
    };

    Vector &Manifold::calA(const Variable &x, const Vector &DLambda, const Problem *prob, Vector *result) const
    {
        return ObtainNorVerIntr(x, DLambda, result);
        /*for debug*/
        /*testCalACalAadj(x, DLambda, ELambda); */
    };

    Vector &Manifold::calAadj(const Variable &x, const Vector &d, const Problem *prob, Vector *result) const
    {
        return ObtainNorVerExtr(x, d, result);
    };
/*
    void RPG::testCalACalAadj(Variable *x, Vector *DLambda, realdp *ELambda)
    {
        Prob->GetDomain()->ObtainNorVerIntr(x, DLambda, ELambda);
        realdp *tmpp = new realdp[dimNorVec];
        for(integer i = 0; i < dimNorVec; i++)
            tmpp[i] = genrandnormal();
        std::cout << "<A d, t>:" << dot_(&dimNorVec, ELambda, &GLOBAL::IONE, tmpp, &GLOBAL::IONE) << std::endl;

        Vector *Vtmp = DLambda->ConstructEmpty();
        calAadj(x, tmpp, Vtmp);
        std::cout << "<d, A^* t>:" << Prob->GetDomain()->Metric(nullptr, Vtmp, DLambda) << std::endl;
        delete Vtmp;
        delete[] tmpp;
    }
*/
}; /*end of ROPTLIB namespace*/
