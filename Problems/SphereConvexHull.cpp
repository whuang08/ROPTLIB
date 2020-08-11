
#include "Problems/SphereConvexHull.h"

/*Define the namespace*/
namespace ROPTLIB{

	SphereConvexHull::SphereConvexHull(const Manifold *inMani, Vector *inW, integer inlengthW, SolversNSMSub *insolver, Vector &(SolversNSMSub::*inHv)(const Vector &v, Vector *result))
	{
		Mani = inMani;
		W = inW;
		lengthW = inlengthW;
		solver = insolver;
		Hv = inHv;
        
        NumGradHess = false;
	};

	SphereConvexHull::~SphereConvexHull(void)
	{
	};

	realdp SphereConvexHull::f(const Variable &x) const
	{
        Vector xsq = x.GetHadamardProduct(x);
        const realdp *xsqPtr = xsq.ObtainReadData();
        
        Vector Wxsq(W[0].Getrow()), PWxsq(W[0].Getrow());
        Wxsq.SetToZeros();
        for(integer i = 0; i < lengthW; i++)
        {
            Wxsq = xsqPtr[i] * W[i] + Wxsq;
        }
        if(Hv == nullptr)
        {
            PWxsq = Wxsq;
        }
        else
        {
            (solver->*Hv)(Wxsq, &PWxsq);
        }
        x.AddToFields("Wxsq", Wxsq);
        x.AddToFields("PWxsq", PWxsq);
        
        return Wxsq.DotProduct(PWxsq);
	};

	Vector &SphereConvexHull::EucGrad(const Variable &x, Vector *result) const
	{
        Vector PWxsq = x.Field("PWxsq");
        
        realdp *resultptr = result->ObtainWriteEntireData();
        const realdp *xptr = x.ObtainReadData();
        for(integer i = 0; i < result->Getlength(); i++)
        {
            resultptr[i] = 4 * xptr[i] * W[i].DotProduct(PWxsq);
        }
        return *result;
	};

	Vector &SphereConvexHull::EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const
	{
        Vector PWxsq = x.Field("PWxsq");
        realdp *resultptr = result->ObtainWriteEntireData();
        const realdp *etaxptr = etax.ObtainReadData();
        for(integer i = 0; i < result->Getlength(); i++)
        {
            resultptr[i] = 4 * etaxptr[i] * W[i].DotProduct(PWxsq);
        }
        Vector xeta = x.GetHadamardProduct(etax);
        const realdp *xetaptr = xeta.ObtainReadData();
        Vector Wxeta(W[0].Getrow()), PWxeta(W[0].Getrow());
        Wxeta.SetToZeros();
        for(integer i = 0; i < lengthW; i++)
        {
            Wxeta = xetaptr[i] * W[i] + Wxeta;
        }
        if(Hv == nullptr)
        {
            PWxeta = Wxeta;
        }
        else
        {
            (solver->*Hv)(Wxeta, &PWxeta);
        }
        
        const realdp *xptr = x.ObtainReadData();
        for(integer i = 0; i < result->Getlength(); i++)
        {
            resultptr[i] += 8 * xptr[i] * W[i].DotProduct(PWxeta);
        }
        
        return *result;
	};

}; /*end of ROPTLIB namespace*/
