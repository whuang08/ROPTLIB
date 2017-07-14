/*
This file defines the class for the Rayleigh Quotient problem on the tangent space of a given manifold
min_{v \in \T_x M} g_x(v, H v), where M is an input manifold. H v is the action of the Hessian of an input problem.
This problem is used to compute the smallest or largest eigenvalue of the action of the Riemannian Hessian of
some problems

Problem --> SphereTxRQ

---- WH
*/

#ifndef SPHERETXRQ_H
#define SPHERETXRQ_H

#include "Manifolds/SphereTx/SphereTx.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Others/def.h"
#include "Others/MyMatrix.h"

/*Define the namespace*/
namespace ROPTLIB{

	class SphereTxRQ : public Problem{
	public:
		SphereTxRQ(Manifold *inmani, Variable *inroot, Problem *inprob, bool inismin = true);
		void SetMinorMax(bool inismin);
		virtual ~SphereTxRQ();
		virtual double f(Variable *x) const;

		virtual void EucGrad(Variable *x, Vector *egf) const;
		virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

		Manifold *mani;
		Problem *prob;
		Variable *root;
		bool ismin;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of GRASSRQ_H
