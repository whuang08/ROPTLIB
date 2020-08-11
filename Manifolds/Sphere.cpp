
#include "Manifolds/Sphere.h"

/*Define the namespace*/
namespace ROPTLIB{

	Sphere::Sphere(integer inn) :Stiefel(inn, 1)
	{
		name.assign("Sphere");
	};

	Sphere::~Sphere(void)
	{
	};

	/* choose qf retraction, parallelization and intrinsic approach and no householder reflections */
	void Sphere::ChooseParamsSet1(void)
	{
		metric = STIE_EUCLIDEAN;
		retraction = STIE_QF;
		VecTran = STIE_PARALLELIZATION;
		IsIntrApproach = true;
		HasHHR = false;
	};
	/* choose exponential map, parallel translation and extrinsic approach and no householder reflections
	 Even though the Householder reflections are not used, the locking condition is satisfied. */
	void Sphere::ChooseParamsSet2(void)
	{
		metric = STIE_EUCLIDEAN;
		retraction = STIE_EXP;
		VecTran = STIE_PARALLELTRANSLATION;
		IsIntrApproach = false;
		HasHHR = false;
	};

	/* choose qf, parallel translation and extrinsic approach and no householder reflections
	The locking conidition is not satisfied */
	void Sphere::ChooseParamsSet3(void)
	{
		metric = STIE_EUCLIDEAN;
		retraction = STIE_QF;
		VecTran = STIE_PARALLELTRANSLATION;
		IsIntrApproach = false;
		HasHHR = false;
	};

	/* choose qf, parallel translation and extrinsic approach and no householder reflections
	Beta \neq 1 is used and the locking conidition is satisfied */
	void Sphere::ChooseParamsSet4(void)
	{
		metric = STIE_EUCLIDEAN;
		retraction = STIE_QF;
		VecTran = STIE_PARALLELTRANSLATION;
		IsIntrApproach = false;
		HasHHR = false;
	};

	Variable &Sphere::ExpRetraction(const Variable &x, const Vector &etax, Variable *result) const
	{
		realdp normetax = sqrt(Metric(x, etax, etax));
		VectorLinearCombination(x, cos(normetax), x, sin(normetax) / normetax, etax, result);
		realdp normresult = sqrt(Metric(x, *result, *result));
		return ScalarTimesVector(x, static_cast<realdp> (1) / normresult, *result, result);
	};

	Vector &Sphere::ExpcoTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
	{
		realdp xiytx = Metric(x, x, xiy);
		realdp xiytetax = Metric(x, xiy, etax);
		realdp normetax = sqrt(Metric(x, etax, etax));
		realdp sinnormetax = sin(normetax);
		realdp cosnormetax = cos(normetax);
		VectorLinearCombination(x, sinnormetax / normetax, xiy, (xiytetax * cosnormetax / normetax
			- xiytx * sinnormetax - xiytetax * sinnormetax / normetax / normetax) / normetax, etax, result);
		return Projection(x, *result, result);
	};

	Vector &Sphere::ExpDiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir) const
	{
		realdp etaxtxix = Metric(x, etax, xix);
		realdp normetax = sqrt(Metric(x, etax, etax));
		realdp sinnormetax = sin(normetax);
        VectorLinearCombination(x, -sinnormetax * etaxtxix / normetax, x, sinnormetax / normetax, xix, result);
        ScalarVectorAddVector(x, (cos(normetax) - sinnormetax / normetax) * etaxtxix / normetax / normetax, etax, *result, result);
        

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
	};

	Vector &Sphere::ExpVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const
	{
        Vector xpy(x); xpy.AlphaXaddThis(1, y); /*xpy = x + y; */
		realdp xynormsq = Metric(x, xpy, xpy);
		realdp xixty = Metric(x, xix, y);
        *result = xix; result->AlphaXaddThis((static_cast<realdp> (-2) * xixty / xynormsq), xpy); /* result = (static_cast<realdp> (-2) * xixty / xynormsq) * xpy + xix; */
        return *result;
	};

	Vector &Sphere::ExpInverseVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
	{
        Vector xpy(x); xpy.AlphaXaddThis(1, y); /*xpy = x + y; */
        realdp xynormsq = Metric(x, xpy, xpy);
		realdp xiytx = Metric(x, xiy, x);
        *result = xiy; result->AlphaXaddThis((static_cast<realdp> (-2) * xiytx / xynormsq), xpy);/*(static_cast<realdp> (-2) * xiytx / xynormsq) * xpy + xiy;*/
        return *result;
	};

	LinearOPE &Sphere::ExpHInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, integer start, integer end, LinearOPE *result) const
	{
        Vector xpy(x); xpy.AlphaXaddThis(1, y); /*xpy = x + y; */
        integer ell = Hx.Getrow();
		integer length = etax.Getlength();
		const realdp *M = Hx.ObtainReadData();
		realdp *Hxpy = new realdp[ell];
		const realdp *xpyTV = xpy.ObtainReadData();

		/* Hxpy <- M(: start : start + length - 1) * xpyTV, details: www.netlib.org/lapack/explore-html/dc/da8/gemv_8f.html */
		gemv_(GLOBAL::N, &ell, &length, &GLOBAL::DONE, const_cast<realdp *> (M + start * ell), &ell, const_cast<realdp *> (xpyTV), &GLOBAL::IONE, &GLOBAL::DZERO, Hxpy, &GLOBAL::IONE);

		realdp scalar = static_cast<realdp> (-2) / Metric(x, xpy, xpy);
        *result = Hx;
		const realdp *xv = x.ObtainReadData();
		realdp *resultL = result->ObtainWritePartialData();
		/* resultL(:, start : start + length - 1) <- scalar * Hxpy * xv^T + resultL(:, start : start + length - 1),
		details: www.netlib.org/lapack/explore-html/dc/da8/ger_8f.html */
		ger_(&length, &ell, &scalar, Hxpy, &GLOBAL::IONE, const_cast<realdp *> (xv), &GLOBAL::IONE, resultL + start * ell, &ell);

		delete[] Hxpy;
        return *result;
	};

	LinearOPE &Sphere::ExpTranH(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, integer start, integer end, LinearOPE *result) const
	{
        Vector xpy(x); xpy.AlphaXaddThis(1, y); /*xpy = x + y; */
        integer ell = Hx.Getrow();
		integer length = etax.Getlength();
		const realdp *M = Hx.ObtainReadData();
		realdp *Hty = new realdp[ell];
		const realdp *yv = y.ObtainReadData();

		/* Hty <- M(start : start + length - 1, :)^T * yv, details: www.netlib.org/lapack/explore-html/dc/da8/gemv_8f.html */
		gemv_(GLOBAL::T, &length, &ell, &GLOBAL::DONE, const_cast<realdp *> (M + start), &ell, const_cast<realdp *> (yv), &GLOBAL::IONE, &GLOBAL::DZERO, Hty, &GLOBAL::IONE);

		realdp scalar = static_cast<realdp> (-2) / Metric(x, xpy, xpy);
		const realdp *xpyTV = xpy.ObtainReadData();
        *result = Hx;
		realdp *resultL = result->ObtainWritePartialData();
		/* resultL(start : start + length - 1, :) <- scalar * xpyTV * Hty^T + resultL(start : start + length - 1, :),
		details: www.netlib.org/lapack/explore-html/dc/da8/ger_8f.html */
		ger_(&length, &ell, &scalar, const_cast<realdp *> (xpyTV), &GLOBAL::IONE, Hty, &GLOBAL::IONE, resultL + start, &ell);

		delete[] Hty;
        return *result;
	};

	LinearOPE &Sphere::ExpTranHInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, LinearOPE *result) const
	{
		ExpHInvTran(x, etax, y, Hx, 0, etax.Getlength(), result);
		return ExpTranH(x, etax, y, *result, 0, etax.Getlength(), result);
	};

	Variable &Sphere::Retraction(const Variable &x, const Vector &etax, Variable *result) const
	{
		if (retraction == STIE_EXP)
		{
			return ExpRetraction(x, etax, result);
		}

		return Stiefel::Retraction(x, etax, result);
	};

	Vector &Sphere::coTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
	{
		if (retraction == STIE_EXP)
		{
			return ExpcoTangentVector(x, etax, y, xiy, result);
		}
		return Stiefel::coTangentVector(x, etax, y, xiy, result);
	};

	Vector &Sphere::DiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir) const
	{
		if (retraction == STIE_EXP)
		{
			return ExpDiffRetraction(x, etax, y, xix, result, IsEtaXiSameDir);
		}
		return Stiefel::DiffRetraction(x, etax, y, xix, result, IsEtaXiSameDir);
	};

	Vector &Sphere::VectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const
	{
		if (VecTran == STIE_PARALLELTRANSLATION)
		{
			return ExpVectorTransport(x, etax, y, xix, result);
		}
		return Stiefel::VectorTransport(x, etax, y, xix, result);
	};

	Vector &Sphere::InverseVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
	{
		if (VecTran == STIE_PARALLELTRANSLATION)
		{
			return ExpInverseVectorTransport(x, etax, y, xiy, result);
		}
		return Stiefel::InverseVectorTransport(x, etax, y, xiy, result);
	};

	LinearOPE &Sphere::HInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, integer start, integer end, LinearOPE *result) const
	{
		if (VecTran == STIE_PARALLELTRANSLATION)
		{
			return ExpHInvTran(x, etax, y, Hx, start, end, result);
		}
		return Stiefel::HInvTran(x, etax, y, Hx, start, end, result);
	};

	LinearOPE &Sphere::TranH(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, integer start, integer end, LinearOPE *result) const
	{
		if (VecTran == STIE_PARALLELTRANSLATION)
		{
			return ExpTranH(x, etax, y, Hx, start, end, result);
		}
		return Stiefel::TranH(x, etax, y, Hx, start, end, result);
	};

	LinearOPE &Sphere::TranHInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, LinearOPE *result) const
	{
		if (VecTran == STIE_PARALLELTRANSLATION)
		{
			return ExpTranHInvTran(x, etax, y, Hx, result);
		}
		return Stiefel::TranHInvTran(x, etax, y, Hx, result);
	};

	void Sphere::SetParams(PARAMSMAP params)
	{
		Stiefel::SetParams(params);
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
