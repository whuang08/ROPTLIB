
#include "Manifolds/SphereTx/SphereTx.h"

/*Define the namespace*/
namespace ROPTLIB{

	SphereTx::SphereTx(Manifold *inmani, Variable *inroot)
	{
		mani = inmani;
		root = inroot;

		// public parameter
		HasHHR = false;

		IsIntrApproach = false;
		UpdBetaAlone = false;

		// Status of locking condition
		HasLockCon = false;

		// strictly speaking, both are extrinsic representations
		// upper one has codimension ExtrinsicDim - IntrinsicDim - 1
		// bottom one has codimension 1
		ExtrinsicDim = mani->GetExtrDim();
		IntrinsicDim = mani->GetIntrDim();
		name.assign("SphereTx");

		if(mani->GetIsIntrinsic())
			EMPTYEXTR = mani->GetEMPTYINTR()->ConstructEmpty();
		else
			EMPTYEXTR = mani->GetEMPTYEXTR()->ConstructEmpty();

		EMPTYINTR = nullptr;
	};

	SphereTx::~SphereTx()
	{
		delete EMPTYEXTR;
	};

	Variable *SphereTx::RandominManifold()
	{
		Variable *result = nullptr;
		if (mani->GetIsIntrinsic())
			result = mani->GetEMPTYINTR()->ConstructEmpty();
		else
			result = mani->GetEMPTYEXTR()->ConstructEmpty();
		result->RandGaussian();
		if (!(mani->GetIsIntrinsic()))
			mani->ExtrProjection(root, result, result);

		double normresult = sqrt(mani->Metric(root, result, result));
		mani->ScaleTimesVector(root, 1.0 / normresult, result, result);
		return result;
	};

	double SphereTx::Metric(Variable *x, Vector *etax, Vector *xix) const
	{
		return mani->Metric(root, etax, xix);
	};

	void SphereTx::Projection(Variable *x, Vector *v, Vector *result) const
	{
		const double *xl = x->ObtainReadData();
		double nume = Metric(x, x, v);
		scalarVectorAddVector(x, -nume, x, v, result);
	};

	void SphereTx::Retraction(Variable *x, Vector *etax, Variable *result, double stepsize) const
	{// exponential mapping
		double norm = sqrt(Metric(x, etax, etax));
		if (norm < std::numeric_limits<double>::epsilon())
			ScaleTimesVector(x, cos(norm), x, result);
		else
			VectorLinearCombination(x, cos(norm), x, sin(norm) / norm, etax, result);

		norm = sqrt(this->Metric(x, result, result));
		this->ScaleTimesVector(x, 1.0 / norm, result, result);
	};

	void SphereTx::coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
	{
		xiy->CopyTo(result);
		printf("The cotangent vector has not been implemented!\n");
	};

	void SphereTx::DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const
	{
		if (IsEtaXiSameDir)
		{
			VectorTransport(x, etax, y, xix, result);

			if (IsEtaXiSameDir && (HasHHR || UpdBetaAlone))
			{
				const double *etaxTV = etax->ObtainReadData();
				const double *xixTV = xix->ObtainReadData();
				double EtatoXi = sqrt(Metric(x, etax, etax) / Metric(x, xix, xix));
				SharedSpace *beta = new SharedSpace(1, 3);
				double *betav = beta->ObtainWriteEntireData();
				betav[0] = sqrt(Metric(x, etax, etax) / Metric(x, result, result)) / EtatoXi;
				betav[1] = Metric(x, etax, etax);
				betav[2] = Metric(x, result, result) * EtatoXi * EtatoXi;
				etax->AddToTempData("beta", beta);

				if (HasHHR)
				{
					Vector *TReta = result->ConstructEmpty();
					result->CopyTo(TReta);
					ScaleTimesVector(x, betav[0] * EtatoXi, TReta, TReta);
					SharedSpace *SharedTReta = new SharedSpace(TReta);
					etax->AddToTempData("betaTReta", SharedTReta);
				}
			}
			return;
		}
		printf("Warning: The differentiated retraction has not been implemented!\n");
		xix->CopyTo(result);
	};

	double SphereTx::Beta(Variable *x, Vector *etax) const
	{
		return 1;
	};

	void SphereTx::VectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const
	{
		if (!etax->TempDataExist("xdydn2"))
		{
			Vector *xdy = x->ConstructEmpty();
			SharedSpace *Sharedxdy = new SharedSpace(xdy);
			VectorAddVector(x, x, y, xdy);
			ScaleTimesVector(x, 1.0 / Metric(x, xdy, xdy), xdy, xdy);
			etax->AddToTempData("xdydn2", Sharedxdy);
		}
		const SharedSpace *Sharedxdydn2 = etax->ObtainReadTempData("xdydn2");
		Vector *xdydn2 = Sharedxdydn2->GetSharedElement();
		scalarVectorAddVector(x, -2.0 * Metric(x, xix, y), xdydn2, xix, result);
	};

	void SphereTx::InverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
	{
		if (!etax->TempDataExist("xdydn2"))
		{
			Vector *xdy = x->ConstructEmpty();
			SharedSpace *Sharedxdy = new SharedSpace(xdy);
			VectorAddVector(x, x, y, xdy);
			ScaleTimesVector(x, 1.0 / Metric(x, xdy, xdy), xdy, xdy);
			etax->AddToTempData("xdydn2", Sharedxdy);
		}
		const SharedSpace *Sharedxdydn2 = etax->ObtainReadTempData("xdydn2");
		Vector *xdydn2 = Sharedxdydn2->GetSharedElement();
		scalarVectorAddVector(x, -2.0 * Metric(x, xiy, x), xdydn2, xiy, result);
	};

	void SphereTx::ObtainIntr(Variable *x, Vector *etax, Vector *result) const
	{
		printf("Routine of obtaining intrinsic representations has not been done!\n");
		etax->CopyTo(result);
	};

	void SphereTx::ObtainExtr(Variable *x, Vector *intretax, Vector *result) const
	{
		printf("Routine of obtaining extrinsic representations has not been done!\n");
		intretax->CopyTo(result);
	};

	void SphereTx::IntrProjection(Variable *x, Vector *v, Vector *result) const
	{
		v->CopyTo(result);
	};

	void SphereTx::ExtrProjection(Variable *x, Vector *v, Vector *result) const
	{
		const double *xl = x->ObtainReadData();
		double nume = Metric(x, x, v);
		scalarVectorAddVector(x, -nume, x, v, result);
	};

	void SphereTx::CheckParams(void) const
	{
		Manifold::CheckParams();
		printf("%s PARAMETERS:\n", name.c_str());
		printf("Tangent space :%15s,\n", mani->GetName().c_str());
	};

	void SphereTx::EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const
	{
		if (prob->GetUseHess())
		{
			Vector *segf = egf->ConstructEmpty();
			segf->NewMemoryOnWrite(); // I don't remember the reason. It seems to be required.
			egf->CopyTo(segf);
			SharedSpace *Sharedegf = new SharedSpace(segf);
			x->AddToTempData("EGrad", Sharedegf);
		}
		ExtrProjection(x, egf, gf);
	};

	void SphereTx::EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const
	{
		const SharedSpace *Sharedegf = x->ObtainReadTempData("EGrad");
		Vector *egfVec = Sharedegf->GetSharedElement();
		double a1 = Metric(x, egfVec, x);
		scalarVectorAddVector(x, -a1, etax, exix, xix);
		ExtrProjection(x, xix, xix);
	};
}; /*end of ROPTLIB namespace*/
