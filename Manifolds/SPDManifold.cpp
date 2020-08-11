
#include "Manifolds/SPDManifold.h"

/*Define the namespace*/
namespace ROPTLIB{

	SPDManifold::SPDManifold(integer inn)
	{
		n = inn;
		IsIntrApproach = true;
		HasHHR = false;
		name.assign("SPDManifold");
		IntrinsicDim = n * (n + 1) / 2;
		ExtrinsicDim = n * n;
		metric = SPDAFFINEINVARIANCE;
        retr = SPDSECONDORDER;
        VecTran = SPDVTPARA;
        EMPTYEXTR = Vector(n, n);
        EMPTYINTR = Vector(IntrinsicDim);
	};

	SPDManifold::~SPDManifold(void)
	{
	};

	void SPDManifold::ChooseParamsSet1(void)
	{
        IsIntrApproach = true;
		metric = SPDAFFINEINVARIANCE;
        retr = SPDSECONDORDER;
        VecTran = SPDVTPARA;
	};

	void SPDManifold::ChooseParamsSet2(void)
	{
        IsIntrApproach = true;
		metric = SPDEUCLIDEAN;
        retr = SPDSECONDORDER;
        VecTran = SPDVTPARA;
	}

    void SPDManifold::ChooseParamsSet3(void)
    {
        IsIntrApproach = false;
        metric = SPDAFFINEINVARIANCE;
        retr = SPDEXP;
        VecTran = SPDPARATRAN;
    };

    void SPDManifold::ChooseParamsSet4(void)
    {
        IsIntrApproach = false;
        metric = SPDEUCLIDEAN;
        retr = SPDEXP;
        VecTran = SPDPARATRAN;
    }

    realdp SPDManifold::Metric(const Variable &x, const Vector &etax, const Vector &xix) const
    {
        if(IsIntrApproach || metric == SPDEUCLIDEAN)
            return etax.DotProduct(xix);
        
        if (metric == SPDAFFINEINVARIANCE)
        {
            Vector tmp = ((x % etax) / x);
            return tmp.DotProduct(xix);
        }
        
        printf("Warning: SPDManifold::metric has not been done!\n");
        return etax.DotProduct(xix);
    };

    Variable SPDManifold::RandominManifold(void) const
    {
        Vector tmp(n, n); tmp.RandGaussian();
        return tmp * tmp.GetTranspose();
    };

	void SPDManifold::CheckParams(void) const
	{
        std::string SPDMetricnames[SPDMETRICLENGTH] = { "EUCLIDEAN", "AFFINV" };
        std::string SPDRetractionnames[SPDRETRACTIONLENGTH] = { "SPDEXP", "SPDSECONDORDER" };
        std::string SPDVecTrannames[SPDVECTORTRANSPORTLENGTH] = { "SPDPARATRAN", "SPDVTPARA" };
		Manifold::CheckParams();
		printf("%s PARAMETERS:\n", name.c_str());
		printf("row           :%15d,\t", n);
		printf("col           :%15d\n", n);
        printf("metric        :%15s,\t", SPDMetricnames[metric].c_str());
        printf("Retraction    :%15s,\n", SPDRetractionnames[retr].c_str());
        printf("VecTran       :%15s,\n", SPDVecTrannames[VecTran].c_str());
	};

	Vector &SPDManifold::ExtrProjection(const Variable &x, const Vector &etax, Vector *result) const
	{
        Vector inetax = etax;
		const realdp *etaxptr = inetax.ObtainReadData();
		realdp *resultptr = result->ObtainWriteEntireData();
		for (integer i = 0; i < n; i++)
		{
			resultptr[i + i * n] = etaxptr[i + i * n];
			for (integer j = i + 1; j < n; j++)
			{
				resultptr[i + j * n] = (etaxptr[i + j * n] + etaxptr[j + i * n]) * static_cast<realdp> (0.5);
				resultptr[j + i * n] = resultptr[i + j * n];
			}
		}
        return *result;
	};

	Vector &SPDManifold::EucGradToGrad(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const
	{
		if (metric == SPDAFFINEINVARIANCE)
			return EucGradToGradAF(x, egf, prob, result);
        
		if (metric == SPDEUCLIDEAN)
			return EucGradToGradEuc(x, egf, prob, result);

        printf("Warning: the Riemannian metric is not valid! Use the affine invariance metric instead!");
        return EucGradToGradAF(x, egf, prob, result);
	};

	Vector &SPDManifold::EucHvToHv(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const
	{
		if (metric == SPDAFFINEINVARIANCE)
			return EucHvToHvAF(x, etax, exix, prob, result);
        
        if (metric == SPDEUCLIDEAN)
            return EucHvToHvEuc(x, etax, exix, prob, result);
        
        printf("Warning: the Riemannian metric is not valid! Use the affine invariance metric instead!");
        return EucHvToHvAF(x, etax, exix, prob, result);
	};

	Vector &SPDManifold::ObtainIntr(const Variable &x, const Vector &etax, Vector *result) const
	{
		if (metric == SPDAFFINEINVARIANCE)
			return ObtainIntrAF(x, etax, result);
        
        if (metric == SPDEUCLIDEAN)
            return ObtainIntrEuc(x, etax, result);

        printf("Warning: the Riemannian metric is not valid! Use the affine invariance metric instead!");
        return ObtainIntrAF(x, etax, result);
	};

	Vector &SPDManifold::ObtainExtr(const Variable &x, const Vector &intretax, Vector *result) const
	{
		if (metric == SPDAFFINEINVARIANCE)
			return ObtainExtrAF(x, intretax, result);
        
        if (metric == SPDEUCLIDEAN)
            return ObtainExtrEuc(x, intretax, result);
        
        printf("Warning: the Riemannian metric is not valid! Use the affine invariance metric instead!");
        return ObtainExtrAF(x, intretax, result);
	};

	realdp SPDManifold::Dist(const Variable &x1, const Variable &x2) const
	{
		if (metric == SPDAFFINEINVARIANCE)
			return DistAF(x1, x2);
        
        if (metric == SPDEUCLIDEAN)
            return DistEuc(x1, x2);
        
        printf("Warning: the Riemannian metric is not valid! Use the affine invariance metric instead!");
        return DistAF(x1, x2);
	};

    Variable &SPDManifold::Retraction(const Variable &x, const Vector &etax, Variable *result) const
    {
        if (retr == SPDEXP)
            return ExpRetraction(x, etax, result);
        
        if (retr == SPDSECONDORDER)
            return SecOrdRetraction(x, etax, result);
        
        printf("Warning: the Riemannian retraction is not valid! Use the second order retraction instead!");
        return SecOrdRetraction(x, etax, result);
    };

    Vector &SPDManifold::DiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir) const
    {
        if (retr == SPDEXP)
            return DiffExpRetraction(x, etax, y, xix, result, IsEtaXiSameDir);
        
        if (retr == SPDSECONDORDER)
            return DiffSecOrdRetraction(x, etax, y, xix, result, IsEtaXiSameDir);
        
        printf("Warning: the Riemannian retraction is not valid! Use the second order retraction instead!");
        return DiffSecOrdRetraction(x, etax, y, xix, result, IsEtaXiSameDir);
    };

    Variable &SPDManifold::ExpRetraction(const Variable &x, const Vector &etax, Variable *result) const
    {
        if (!x.FieldsExist("_L"))
        {
            x.CholDecom();
        }
        Vector exetax(EMPTYEXTR);
        if(IsIntrApproach)
            ObtainExtr(x, etax, &exetax);
        else
            exetax = etax;
        
        Vector LiE = exetax.TriangleLinSol(x.Field("_L"));
        Vector LiELiT = LiE.GetTranspose().TriangleLinSol(x.Field("_L"));
        Vector tmp(n, n);
        tmp.AlphaABaddBetaThis(1, x.Field("_L"), GLOBAL::N, LiELiT.ExpSym(), GLOBAL::N, 0);
        result->AlphaABaddBetaThis(1, tmp, GLOBAL::N, x.Field("_L"), GLOBAL::T, 0);
        realdp *resultptr = result->ObtainWritePartialData();
        
        realdp tmpv = 0;
        for(integer i = 0; i < n; i++)
        {
            for(integer j = i + 1; j < n; j++)
            {
                tmpv = (resultptr[i + j * n] + resultptr[j + i * n]) / 2;
                resultptr[i + j * n] = tmpv;
                resultptr[j + i * n] = tmpv;
            }
        }
        
        if (!result->FieldsExist("_L"))
        {
            result->CholDecom();
        }
        return *result;
    };

    Vector &SPDManifold::DiffExpRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir) const
    {
        printf("Warning: SPDManifold::DiffExpRetraction has not been implemented yet!\n");
        return Manifold::DiffRetraction(x, etax, y, xix, result, IsEtaXiSameDir);
    };

	Variable &SPDManifold::SecOrdRetraction(const Variable &x, const Vector &etax, Variable *result) const
	{ /* result = x + etx + 0.5 * etax * x^{-1} * etax */
        if (!x.FieldsExist("_L"))
        {
            x.CholDecom();
        }
        Vector exetax(EMPTYEXTR);
        if(IsIntrApproach)
            ObtainExtr(x, etax, &exetax);
        else
            exetax = etax;
        
        Vector LiE = exetax.TriangleLinSol(x.Field("_L"));
        
        result->AlphaABaddBetaThis(0.5, LiE, GLOBAL::T, LiE, GLOBAL::N, 0);
        
        result->AlphaXaddThis(1, exetax);
        result->AlphaXaddThis(1, x);
        
        if (!result->FieldsExist("_L"))
        {
            result->CholDecom();
        }
        return *result;
	};

	Vector &SPDManifold::DiffSecOrdRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir) const
	{ /* xix + 0.5 xix x^{-1} etax + 0.5 etax x^{-1} xix */
        realdp nxix = std::sqrt(Metric(x, xix, xix));

        if (!x.FieldsExist("_L"))
        {
            x.CholDecom();
        }
        Vector exetax(EMPTYEXTR), exxix(EMPTYEXTR);
        if(IsIntrApproach)
        {
            ObtainExtr(x, etax, &exetax);
            ObtainExtr(x, xix, &exxix);
        }
        else
        {
            exetax = etax;
            exxix = xix;
        }
        Vector LiE = exetax.TriangleLinSol(x.Field("_L"));
        Vector LiX = exxix.TriangleLinSol(x.Field("_L"));
        Vector exresult(EMPTYEXTR);
        exresult.AlphaABaddBetaThis(1, LiE, GLOBAL::T, LiX, GLOBAL::N, 0);
        realdp *exresultptr = exresult.ObtainWritePartialData();
        for (integer i = 0; i < n; i++)
        {
            for (integer j = i + 1; j < n; j++)
            {
                exresultptr[j + i * n] = (exresultptr[j + i * n] + exresultptr[i + j * n]) * static_cast<realdp> (0.5);
                exresultptr[i + j * n] = exresultptr[j + i * n];
            }
        }
        exresult.AlphaXaddThis(1, exxix);
        
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

	Vector &SPDManifold::coTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
	{
		printf("SPDManifold::coTangentVector has not been done!\n");
        return Manifold::coTangentVector(x, etax, y, xiy, result);
	};

	realdp SPDManifold::Beta(const Variable &x, const Vector &etax) const
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

	Vector &SPDManifold::EucGradToGradAF(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const
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
        Vector tmp(n, n);
        tmp.AlphaABaddBetaThis(1, x, GLOBAL::T, egf, GLOBAL::N, 0);
        result->AlphaABaddBetaThis(1, tmp, GLOBAL::N, x, GLOBAL::N, 0);
        return *result;
	};

	Vector &SPDManifold::EucHvToHvAF(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const
	{
        Vector tmp = x * x.Field("EGrad") * etax;
        *result = x * exix * x + (tmp + tmp.GetTranspose()) / 2;
        return *result;
	};

	Vector &SPDManifold::ObtainIntrAF(const Variable &x, const Vector &etax, Vector *result) const
	{/* L^{-1} etax L^{-T}, where x = L L^T, etax is assumed to be a symmetric matrix. */
        
        if (!x.FieldsExist("_L"))
        {
            x.CholDecom();
        }
        
        Vector LiE = etax.TriangleLinSol(x.Field("_L"));
        Vector LiELiT = LiE.GetTranspose().TriangleLinSol(x.Field("_L"));
        const realdp *LiELiTptr = LiELiT.ObtainReadData();
        
		realdp *resultptr = result->ObtainWriteEntireData();
		integer idx = 0;
		realdp r2 = static_cast<realdp> (sqrt(2.0));
		for (integer i = 0; i < n; i++)
		{
			resultptr[idx] = LiELiTptr[i + i * n];
			idx++;
		}

		for (integer i = 0; i < n; i++)
		{
			for (integer j = i + 1; j < n; j++)
			{
				resultptr[idx] = LiELiTptr[j + i * n] * r2;
				idx++;
			}
		}
        return *result;
	};

	Vector &SPDManifold::ObtainExtrAF(const Variable &x, const Vector &intretax, Vector *result) const
	{
        if (!x.FieldsExist("_L"))
        {
            x.CholDecom();
        }
        
		const realdp *intretaxptr = intretax.ObtainReadData();
		realdp *resultptr = result->ObtainWriteEntireData();

		integer idx = 0;
		realdp r2 = static_cast<realdp> (sqrt(2.0));
		for (integer i = 0; i < n; i++)
		{
			resultptr[i + i * n] = intretaxptr[idx];
			idx++;
		}

		for (integer i = 0; i < n; i++)
		{
			for (integer j = i + 1; j < n; j++)
			{
				resultptr[j + i * n] = intretaxptr[idx] / r2;
				resultptr[i + j * n] = resultptr[j + i * n];
				idx++;
			}
		}
        
        Vector tmp(n, n); tmp.AlphaABaddBetaThis(1, x.Field("_L"), GLOBAL::N, *result, GLOBAL::N, 0);
        result->AlphaABaddBetaThis(1, tmp, GLOBAL::N, x.Field("_L"), GLOBAL::T, 0);
        
        return *result;
	};

	realdp SPDManifold::DistAF(const Variable &x1, const Variable &x2) const
	{ /* dis(x, y) = |log(Lx^{-1} y Lx^{-T})|_F */
		
        if (!x1.FieldsExist("_L"))
        {
            x1.CholDecom();
        }
        
        if (!x2.FieldsExist("_L"))
        {
            x2.CholDecom();
        }
        Vector tmp = x2.Field("_L").TriangleLinSol(x1.Field("_L"));
        Vector tmp2(EMPTYEXTR); tmp2.AlphaABaddBetaThis(1, tmp, GLOBAL::N, tmp, GLOBAL::T, 0);
        tmp2.EigenDecomSym();
        Vector EigVals = tmp2.Field("_EigVal");
        return EigVals.Fnorm();
	};

	Vector &SPDManifold::EucGradToGradEuc(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const
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
        *result = egf;
        return *result;
	};

	Vector &SPDManifold::EucHvToHvEuc(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const
	{
        *result = exix;
        return *result;
	};

    Vector &SPDManifold::VectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const
    {
        if (VecTran == SPDPARATRAN)
            return ParallelTranslation(x, etax, y, xix, result);
        
        if (VecTran == SPDVTPARA)
            return VectorTransportParallelization(x, etax, y, xix, result);

        printf("Warning: the VectorTransport is not valid! Use vector transport by parallelization instead!");
        return VectorTransportParallelization(x, etax, y, xix, result);
    };

    Vector &SPDManifold::ParallelTranslation(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const
    {
        if (!x.FieldsExist("_L"))
        {
            x.CholDecom();
        }
        Vector exetax(EMPTYEXTR), exxix(EMPTYEXTR);
        if(IsIntrApproach)
        {
            ObtainExtr(x, etax, &exetax);
            ObtainExtr(x, xix, &exxix);
        }
        else
        {
            exetax = etax;
            exxix = xix;
        }
        
        Vector LiE = exetax.TriangleLinSol(x.Field("_L"));
        Vector LiELiThalf = LiE.GetTranspose().TriangleLinSol(x.Field("_L"));
        LiELiThalf.ScalarTimesThis(0.5);
        Vector tmp = x.Field("_L") * LiELiThalf.ExpSym().TriangleLinSol(x.Field("_L"), GLOBAL::T).GetTranspose();
        LiE.AlphaABaddBetaThis(0, tmp, GLOBAL::N, exxix, GLOBAL::N, 0); /*LiE is for temporary storage*/
        result->AlphaABaddBetaThis(1, LiE, GLOBAL::N, tmp, GLOBAL::T, 0); /* *result = tmp * exxix * tmp.GetTranspose(); */
        return *result;
    };
        
    Vector &SPDManifold::VectorTransportParallelization(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const
    {
        if(IsIntrApproach)
        {
            *result = xix;
            return *result;
        }
        Vector tmp(EMPTYEXTR); ObtainIntr(x, xix, &tmp);
        return ObtainExtr(y, tmp, result);
    };

	Vector &SPDManifold::ObtainIntrEuc(const Variable &x, const Vector &etax, Vector *result) const
	{
		const realdp *etaxTV = etax.ObtainReadData();
		realdp *resultTV = result->ObtainWriteEntireData();
		integer idx = 0;
		realdp r2 = static_cast<realdp> (sqrt(2.0));
		for (integer i = 0; i < n; i++)
		{
			resultTV[idx] = etaxTV[i + i * n];
			idx++;
		}

		for (integer i = 0; i < n; i++)
		{
			for (integer j = i + 1; j < n; j++)
			{
				resultTV[idx] = etaxTV[j + i * n] * r2;
				idx++;
			}
		}
        return *result;
	};

	Vector &SPDManifold::ObtainExtrEuc(const Variable &x, const Vector &intretax, Vector *result) const
	{
		const realdp *intretaxTV = intretax.ObtainReadData();
		realdp *resultTV = result->ObtainWriteEntireData();

		integer idx = 0;
		realdp r2 = static_cast<realdp> (sqrt(2.0));
		for (integer i = 0; i < n; i++)
		{
			resultTV[i + i * n] = intretaxTV[idx];
			idx++;
		}

		for (integer i = 0; i < n; i++)
		{
			for (integer j = i + 1; j < n; j++)
			{
				resultTV[j + i * n] = intretaxTV[idx] / r2;
				resultTV[i + j * n] = resultTV[j + i * n];
				idx++;
			}
		}
        return *result;
	};

	realdp SPDManifold::DistEuc(const Variable &x1, const Variable &x2) const
	{
        return (x1 - x2).Fnorm();
	};
}; /*end of ROPTLIB namespace*/
