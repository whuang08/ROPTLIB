
#include "Problems/SphereTxRQ/SphereTxRQ.h"

/*Define the namespace*/
namespace ROPTLIB{

	SphereTxRQ::SphereTxRQ(Manifold *inmani, Variable *inroot, Problem *inprob, bool inismin)
	{
		mani = inmani;
		prob = inprob;
		root = inroot;
		ismin = inismin;
		prob->SetUseGrad(true);
		prob->SetUseHess(true);
		prob->f(root);
		Vector *gf = nullptr;
		if (mani->GetIsIntrinsic())
			gf = mani->GetEMPTYINTR()->ConstructEmpty();
		else
			gf = mani->GetEMPTYEXTR()->ConstructEmpty();
		prob->Grad(root, gf);
		delete gf;
	};

	void SphereTxRQ::SetMinorMax(bool inismin)
	{
		ismin = inismin;
	};

	SphereTxRQ::~SphereTxRQ(void)
	{
	};

	double SphereTxRQ::f(Variable *x) const
	{
		Variable *Hx = x->ConstructEmpty();
		SharedSpace *SharedHx = new SharedSpace(Hx);
		prob->HessianEta(root, x, Hx);
		double result = mani->Metric(root, x, Hx);
		x->AddToTempData("Hx", SharedHx);
		result = ismin ? result : -result;
		return result;
	};

	void SphereTxRQ::EucGrad(Variable *x, Vector *egf) const
	{
		const SharedSpace *SharedHx = x->ObtainReadTempData("Hx");
		Vector *Hx = SharedHx->GetSharedElement();
		if(ismin)
			Domain->ScaleTimesVector(x, 2.0, Hx, egf);
		else
			Domain->ScaleTimesVector(x, -2.0, Hx, egf);
	};

	void SphereTxRQ::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const
	{
		prob->HessianEta(root, etax, exix);
		if (ismin)
			Domain->ScaleTimesVector(x, 2.0, exix, exix);
		else
			Domain->ScaleTimesVector(x, -2.0, exix, exix);
	};
}; /*end of ROPTLIB namespace*/
