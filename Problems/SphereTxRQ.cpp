
#include "Problems/SphereTxRQ.h"

/*Define the namespace*/
namespace ROPTLIB{

    Vector MinMaxEigValHessian(Variable *X, Manifold *Domain, const Problem *Prob)
    {
        Vector result(2);
        realdp *resultptr = result.ObtainWriteEntireData();
        {
            SphereTx DomainPH(Domain, X);
            SphereTxRQ ProbHess(Domain, X, Prob, true);
            ProbHess.SetDomain(&DomainPH);
            Variable TV0 = DomainPH.RandominManifold();
            RTRNewton RTRNewtonsolver(&ProbHess, &TV0);
            RTRNewtonsolver.Verbose = NOOUTPUT;
            RTRNewtonsolver.Run();
            if(RTRNewtonsolver.Getnormgfgf0() > 1e-6)
            {
                printf("Warning: Rayleigh quotient for computing the smallest eigenvalue of the Hessian is not solved accurately!\n");
            }
            resultptr[0] = RTRNewtonsolver.Getfinalfun();
        }
        {
            SphereTx DomainPH(Domain, X);
            SphereTxRQ ProbHess(Domain, X, Prob, false);
            ProbHess.SetDomain(&DomainPH);
            Variable TV0 = DomainPH.RandominManifold();
            RTRNewton RTRNewtonsolver(&ProbHess, &TV0);
            RTRNewtonsolver.Verbose = NOOUTPUT;
            RTRNewtonsolver.Run();
            if(RTRNewtonsolver.Getnormgfgf0() > 1e-6)
            {
                printf("Warning: Rayleigh quotient for computing the largest eigenvalue of the Hessian is not solved accurately!\n");
            }
            resultptr[1] = - RTRNewtonsolver.Getfinalfun();
        }
        
        return result;
    };

	SphereTxRQ::SphereTxRQ(Manifold *inmani, Variable *inroot, const Problem *inprob, bool inismin)
	{
		mani = inmani;
		prob = inprob;
		root = *inroot;
		ismin = inismin;
		prob->SetUseGrad(true);
		prob->SetUseHess(true);
		prob->f(root);
        
        if(inmani->GetIsIntrinsic())
            TmpTV = inmani->GetEMPTYINTR();
        else
            TmpTV = inmani->GetEMPTYEXTR();
        
        prob->Grad(root, &TmpTV); /*for creating "EGrad" data in root*/
        
        NumGradHess = false;
	};

	void SphereTxRQ::SetMinOrMax(bool inismin)
	{
		ismin = inismin;
	};

	SphereTxRQ::~SphereTxRQ(void)
	{
	};

	realdp SphereTxRQ::f(const Variable &x) const
	{
        prob->HessianEta(root, x, &TmpTV);
        realdp result = mani->Metric(root, TmpTV, x);
        x.AddToFields("Hx", TmpTV);
        result = ismin ? result : - result;
        return result;
	};

	Vector &SphereTxRQ::EucGrad(const Variable &x, Vector *result) const
	{
        *result = x.Field("Hx");
        if(ismin)
        {
            mani->ScalarTimesVector(root, 2, *result, result);
            return *result;
        }
        mani->ScalarTimesVector(root, -2, *result, result);
        
        return *result;
	};

	Vector &SphereTxRQ::EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const
	{
        prob->HessianEta(root, etax, result);
        
        if(ismin)
            return mani->ScalarTimesVector(root, 2, *result, result);
        
        return mani->ScalarTimesVector(root, -2, *result, result);
	};
}; /*end of ROPTLIB namespace*/
