/*
This file defines the class for the Rayleigh Quotient problem on the Grassmann manifold
min_X tr(X^T B X), where B is a symmetric matrix, and X \in Gr(p, n).

Problem --> GrassRQ

---- WH
*/

#ifndef GRASSRQ_H
#define GRASSRQ_H

#include "Manifolds/Grassmann.h"
#include "Problems/Problem.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class GrassRQ : public Problem{
	public:
		GrassRQ(Vector inB, integer inn, integer inp);
		virtual ~GrassRQ();
		virtual realdp f(const Variable &x) const;

		virtual Vector &EucGrad(const Variable &x, Vector *result) const;
		virtual Vector &EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const;

		Vector B;
		integer n;
		integer p;
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of GRASSRQ_H */
