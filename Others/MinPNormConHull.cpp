
#include "Others/MinPNormConHull.h"

/*Define the namespace*/
namespace ROPTLIB{

	realdp MinPNormConHull(const Manifold *Mani, Variable x, Vector *Ys, integer LYs, Vector &outSoln)
	{
		return MinPNormConHull(Mani, x, Ys, LYs, nullptr, nullptr, outSoln);
	};

	/*Compute min_{y in convex hull of gfs, and gfs are tangent vectors at the tangent space at x} ||y|| */
	realdp MinPNormConHull(const Manifold *Mani, Variable x, Vector *Ys, integer LYs, SolversNSMSub *solver, Vector &(SolversNSMSub::*Hv)(const Vector &v, Vector *result), Vector &outSoln)
	{
		/*Use Riemannian method*/
		return MinPNormConHullRMethod(Mani, x, Ys, LYs, solver, Hv, outSoln);
	};

	realdp MinPNormConHullRMethod(const Manifold *Mani, Variable x, Vector *Ys, integer LYs, SolversNSMSub *solver, Vector &(SolversNSMSub::*Hv)(const Vector &v, Vector *result), Vector &outSoln)
	{
		SphereConvexHull subprob(Mani, Ys, LYs, solver, Hv);
		Sphere Domain(LYs);
		subprob.SetDomain(&Domain);
        Variable InitX = Domain.RandominManifold();

		RTRNewton *RTRNewtonsolver = new RTRNewton(&subprob, &InitX);
		RTRNewtonsolver->Stop_Criterion = SM_GRAD_F;
        RTRNewtonsolver->Verbose = NOOUTPUT;
		RTRNewtonsolver->Max_Iteration = 100;
		RTRNewtonsolver->Tolerance = static_cast<realdp> (1e-10);
		RTRNewtonsolver->Run();
		realdp fopt = RTRNewtonsolver->Getfinalfun();
        Vector xopt = RTRNewtonsolver->GetXopt();
        
        outSoln = xopt.Field("Wxsq");
		delete RTRNewtonsolver;
		return fopt;
	};
}; /*end of ROPTLIB namespace*/
