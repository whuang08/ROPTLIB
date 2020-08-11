
#include "Manifolds/Grassmann.h"

/*Define the namespace*/
namespace ROPTLIB{

	Grassmann::Grassmann(integer inn, integer inp)
	{
		HasHHR = false;

		IsIntrApproach = true;

		n = inn;
		p = inp;
		ExtrinsicDim = n * p;
		IntrinsicDim = (n - p) * p;
		name.assign("Grassmann");
        EMPTYEXTR = Vector (n, p);
        EMPTYINTR = Vector (IntrinsicDim);
	};

	Grassmann::~Grassmann(void)
	{
	};

    Variable Grassmann::RandominManifold(void) const
    {
        Variable result(n, p);
        result.RandGaussian();
        result.QRDecom();
        return result.Field("_Q");
    };

	void Grassmann::CheckParams(void) const
	{
		Manifold::CheckParams();
		printf("%s PARAMETERS:\n", name.c_str());
		printf("n             :%15d,\t", n);
		printf("p             :%15d\n", p);
	};

	Vector &Grassmann::ExtrProjection(const Variable &x, const Vector &etax, Variable *result) const
	{
        Vector xTetax(p, p);
        xTetax.AlphaABaddBetaThis(1, x, GLOBAL::T, etax, GLOBAL::N, 0); /* xTetax = x^T * etax*/
        *result = etax;
        result->AlphaABaddBetaThis(-1, x, GLOBAL::N, xTetax, GLOBAL::N, 1); /* result = etax - x * x^T * etax */
        return *result;
	};

	Variable &Grassmann::Retraction(const Variable &x, const Vector &etax, Variable *result) const
	{ /* qf retraction: result = R_x(etax) = qf(x + etax) */
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

	Vector &Grassmann::DiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir) const
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

	Vector &Grassmann::coTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
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

	realdp Grassmann::Beta(const Variable &x, const Vector &etax) const
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

	Vector &Grassmann::EucGradToGrad(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const
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
		return ExtrProjection(x, egf, result);
	};

	Vector &Grassmann::EucHvToHv(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const
	{
        Vector EGrad = x.Field("EGrad");
        Vector tmp(p, p); tmp.AlphaABaddBetaThis(1, x, GLOBAL::T, EGrad, GLOBAL::N, 0); /* tmp = x.GetTranspose() * EGrad; */
        *result = exix; result->AlphaABaddBetaThis(-1, etax, GLOBAL::N, tmp, GLOBAL::N, 1);
        return ExtrProjection(x, *result, result);
	};

	Vector &Grassmann::ObtainIntr(const Variable &x, const Vector &etax, Vector *result) const
	{
        if(!x.FieldsExist("_HHR"))
        {
            x.HHRDecom();
        }
        Vector tmp = etax.HHRMtp(x.Field("_HHR"), x.Field("_tau"), GLOBAL::T, GLOBAL::L);
        
        const realdp *tmpptr = tmp.ObtainReadData();
        realdp *resultptr = result->ObtainWriteEntireData();
        
        for (integer i = 0; i < p; i++)
        {
            integer nmp = n - p;
            copy_(&nmp, const_cast<realdp *>(tmpptr + p + n* i), &GLOBAL::IONE, resultptr + nmp * i, &GLOBAL::IONE);
        }
        
        return *result;
	};

	Vector &Grassmann::ObtainExtr(const Variable &x, const Vector &intretax, Vector *result) const
	{
        if(!x.FieldsExist("_HHR"))
        {
            x.HHRDecom();
        }
        realdp *resultptr = result->ObtainWriteEntireData();
        const realdp *intretaxptr = intretax.ObtainReadData();

        for (integer i = 0; i < p; i++)
        {
            for(integer j = 0; j < p; j++)
            {
                resultptr[j + i * n] = 0;
            }
            integer nmp = n - p;
            copy_(&nmp, const_cast<realdp *> (intretaxptr)+nmp * i, &GLOBAL::IONE, resultptr + p + n* i, &GLOBAL::IONE);
        }

        (*result) = result->HHRMtp(x.Field("_HHR"), x.Field("_tau"), GLOBAL::N, GLOBAL::L);
        
        return *result;
	};
}; /*end of ROPTLIB namespace*/
