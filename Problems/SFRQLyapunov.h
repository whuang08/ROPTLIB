/*
This file defines a problem for finding a low rank approximation for the Lyapunov equation:
A X M + M X A = C C^T,
where A, M, and C are symmetric matrices.
The optimization problem is:
min_{Y \in R_*^{n \times p} / O_p} \tr(Y Y^T A Y Y^T M) - tr(Y Y^T C C^T).

Problem --> SFRQLyapunov

---- WH
*/

#ifndef SFRQLYAPUNOV_H
#define SFRQLYAPUNOV_H

#include "Problems/Problem.h"
#include "Others/def.h"
#include "Manifolds/Element.h"
#include "Manifolds/SymFixedRankQ.h"

/*Define the namespace*/
namespace ROPTLIB{

	class SFRQLyapunov : public Problem{
	public:
		/*define A, M and C in the cost function \tr(Y Y^T A Y Y^T M) - tr(Y Y^T C)*/
        SFRQLyapunov(SparseMatrix &inA, SparseMatrix &inM, Vector inC, integer inp);

		virtual ~SFRQLyapunov();
        
		virtual realdp f(const Variable &x) const;

		virtual Vector &EucGrad(const Variable &x, Vector *result) const;

		virtual Vector &EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const;

        virtual Vector &PreConditioner(const Variable &x, const Vector &eta, Vector *result) const;

		Vector &ActionEH(const Variable &x, const Vector &intreta, Vector *result) const;

		Vector &LinearCG(const Variable x, const Vector &intreta, Vector *result) const;

        
        SparseMatrix *Aptr;
        SparseMatrix *Mptr;
        
		Vector C;

        integer Cp;
		integer n;
		integer p;
	};
}; /*end of ROPTLIB namespace*/

#endif /*SFRQLYAPUNOV_H*/
