/*
This file defines a problem for finding a low rank approximation for the Lyapunov equation:
A X M + M X A = C,
where A, M, and C are symmetric matrices.
The optimization problem is:
min_{Y \in R_*^{n \times p} / O_p} \tr(Y Y^T A Y Y^T M) - tr(Y Y^T C).

Problem --> NSOLyapunov

---- WH
*/

#ifndef NSOLYAPUNOV_H
#define NSOLYAPUNOV_H

#include "Manifolds/NStQOrth/NStQOrth.h"
#include "Manifolds/NStQOrth/NSOVariable.h"
#include "Manifolds/NStQOrth/NSOVector.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class NSOLyapunov : public Problem{
	public:
		/*define A, M and C in the cost function \tr(Y Y^T A Y Y^T M) - tr(Y Y^T C)*/
		NSOLyapunov(double *inA, double *inM, double *inC, integer inn, integer inp);

		virtual ~NSOLyapunov();
		virtual double f(Variable *x) const;

		virtual void EucGrad(Variable *x, Vector *gf) const;

		virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

		double *A;
		double *M;
		double *C;

		integer n;
		integer p;
	};
}; /*end of ROPTLIB namespace*/

#endif NSOLYAPUNOV_H
