
#include "Manifolds/SphereTx.h"

/*Define the namespace*/
namespace ROPTLIB{

	SphereTx::SphereTx(Manifold *inmani, Variable *inroot)
	{
		mani = inmani;
		root = *inroot;

		HasHHR = false;

		IsIntrApproach = false;

		/* strictly speaking, both are extrinsic representations
		 upper one has codimension ExtrinsicDim - IntrinsicDim - 1
		 bottom one has codimension 1 */
		ExtrinsicDim = mani->GetExtrDim();
		IntrinsicDim = mani->GetIntrDim();
		name.assign("SphereTx");
        
        if(mani->GetIsIntrinsic())
            EMPTYEXTR = mani->GetEMPTYINTR();
        else
            EMPTYEXTR = mani->GetEMPTYEXTR();
        
        EMPTYINTR = Vector (0);
	};

	SphereTx::~SphereTx()
	{
	};

	Variable SphereTx::RandominManifold() const
	{
        Vector result(EMPTYEXTR);
		result.RandGaussian();
        mani->Projection(root, result, &result);
		realdp normresult = sqrt(mani->Metric(root, result, result));
        mani->ScalarTimesVector(root, static_cast<realdp> (1) / normresult, result, &result);
		return result;
	};

	realdp SphereTx::Metric(const Variable &x, const Vector &etax, const Vector &xix) const
	{
		return mani->Metric(root, etax, xix);
	};

	Vector &SphereTx::Projection(const Variable &x, const Vector &etax, Vector *result) const
	{
        return ExtrProjection(x, etax, result);
	};

    Vector &SphereTx::ExtrProjection(const Variable &x, const Vector &etax, Vector *result) const
    {
        mani->Projection(root, etax, result);
        realdp nume = mani->Metric(root, x, *result);
        return mani->ScalarVectorAddVector(root, -nume, x, *result, result);
    };

    Vector &SphereTx::ScalarTimesVector(const Variable &x, const realdp &scalar, const Vector &etax, Vector *result) const
    {
        return mani->ScalarTimesVector(root, scalar, etax, result);
    };
    
    Vector &SphereTx::ScalarVectorAddVector(const Variable &x, const realdp &scalar, const Vector &etax, const Vector &xix, Vector *result) const
    {
        return mani->ScalarVectorAddVector(root, scalar, etax, xix, result);
    };

    Vector &SphereTx::VectorLinearCombination(const Variable &x, realdp scalar1, const Vector &etax, realdp scalar2, const Vector &xix, Vector *result) const
    {
        return mani->VectorLinearCombination(root, scalar1, etax, scalar2, xix, result);
    };

	Variable &SphereTx::Retraction(const Variable &x, const Vector &etax, Variable *result) const
	{/* exponential mapping */
        
		realdp norm = sqrt(Metric(x, etax, etax));
		if (norm < std::numeric_limits<realdp>::epsilon())
			mani->ScalarTimesVector(x, cos(norm), x, result);
		else
			mani->VectorLinearCombination(root, cos(norm), x, sin(norm) / norm, etax, result);

		norm = sqrt(this->Metric(x, *result, *result));
        return mani->ScalarTimesVector(root, static_cast<realdp> (1) / norm, *result, result);
	};

	Vector &SphereTx::coTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
	{
		printf("Warning: SphereTx::coTangentVector has not been implemented!\n");
        return Manifold::coTangentVector(x, etax, y, xiy, result);
	};

	Vector &SphereTx::DiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir) const
	{
		if (IsEtaXiSameDir)
		{
            VectorTransport(x, etax, y, xix, result);
            
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
        printf("Warning: SphereTx::DiffRetraction has not been implemented!\n");
        return Manifold::DiffRetraction(x, etax, y, xix, result, IsEtaXiSameDir);
	};

	realdp SphereTx::Beta(const Variable &x, const Vector &etax) const
	{
		return 1;
	};

	Vector &SphereTx::VectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const
	{
        Vector xdy(x); mani->ScalarVectorAddVector(root, 1, y, x, &xdy);
        realdp tmpv = static_cast<realdp> (-2) * mani->Metric(root, xix, y) / mani->Metric(root, xdy, xdy);
        return mani->ScalarVectorAddVector(root, tmpv, xdy, xix, result);
	};

	Vector &SphereTx::InverseVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const
	{
        Vector xdy(y); mani->ScalarVectorAddVector(root, 1, x, y, &xdy);
        realdp tmpv = static_cast<realdp> (-2) * mani->Metric(root, xiy, x) / mani->Metric(root, xdy, xdy);
        return mani->ScalarVectorAddVector(root, tmpv, xdy, xiy, result);
	};

	void SphereTx::CheckParams(void) const
	{
		Manifold::CheckParams();
		printf("%s PARAMETERS:\n", name.c_str());
		printf("Tangent space :%15s,\n", mani->GetName().c_str());
	};

	Vector &SphereTx::EucGradToGrad(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const
	{
        if (prob->GetUseHess())
        {
            x.AddToFields("EGrad", egf);
        }
        return Projection(x, egf, result);
	};

	Vector &SphereTx::EucHvToHv(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const
	{
        Vector EGrad = x.Field("EGrad");
        
        mani->VectorLinearCombination(root, 1, exix, - mani->Metric(root, EGrad, x), etax, result);
        return Projection(x, *result, result);
	};
}; /*end of ROPTLIB namespace*/
