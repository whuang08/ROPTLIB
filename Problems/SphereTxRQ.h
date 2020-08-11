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

#include "Manifolds/SphereTx.h"
#include "Problems/Problem.h"
#include "Solvers/RTRNewton.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

    extern Vector MinMaxEigValHessian(Variable *X, Manifold *Domain, const Problem *Prob);

	class SphereTxRQ : public Problem{
	public:
		SphereTxRQ(Manifold *inmani, Variable *inroot, const Problem *inprob, bool inismin = true);
		void SetMinOrMax(bool inismin);
		virtual ~SphereTxRQ();
		virtual realdp f(const Variable &x) const;

		virtual Vector &EucGrad(const Variable &x, Vector *result) const;
		virtual Vector &EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const;

		Manifold *mani;
		const Problem *prob;
		mutable Variable root;
		bool ismin;
        mutable Vector TmpTV; /*a tangent vector for temporary storage*/
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of SPHERETXRQ_H */
