/*
This file defines the class for 
min_{X \in S^{n-1}} tr((X.^2)^T W^T P W (X.^2)), where P is a N by N symmetric positive definite matrix
and W is a N by n matrix.

Problem --> SphereConvexHull

---- WH
*/

#ifndef SPHERECONVEXHULL_H
#define SPHERECONVEXHULL_H

#include "Manifolds/Sphere.h"
#include "Problems/Problem.h"
#include "Solvers/SolversNSMSub.h"
#include "Solvers/RTRNewton.h"
#include "Solvers/LRBFGS.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

    class SolversLSLPSub;

	class SphereConvexHull : public Problem{
	public:
		/*W is an array of vectors, which are tangent vectors in the tangent space at x, lengthW is the length of W.
		HvRBFGSSub defines a linear mapping: P: v --> Pv. The manifold, Mani, of x and W is not necessary the same as the domain of this problem.*/
		SphereConvexHull(const Manifold *Mani, Vector *W, integer lengthW, SolversNSMSub *solver, Vector &(SolversNSMSub::*Hv)(const Vector &v, Vector *result));
		virtual ~SphereConvexHull();
		virtual realdp f(const Variable &x) const;

		virtual Vector &EucGrad(const Variable &x, Vector *result) const;
		virtual Vector &EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const;

		const Manifold *Mani;
		Vector *W;
		integer lengthW;
		SolversNSMSub *solver;
        Vector &(SolversNSMSub::*Hv)(const Vector &v, Vector *result);
	};

}; /*end of ROPTLIB namespace*/
#endif /* end of GRASSRQ_H */
