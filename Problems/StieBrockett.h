/*
This file defines the class for the problem
min_X tr(X^T B X D), where B is a symmetric matrix, D is a diagonal matrix and X \in St(p, n).

Problem --> StieBrockett

---- WH
*/

#ifndef STIEBROCKETT_H
#define STIEBROCKETT_H

#include "Manifolds/Stiefel.h"
#include "Problems/Problem.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class StieBrockett : public Problem{
	public:
		StieBrockett(Vector inB, Vector inD);
		virtual ~StieBrockett();
		virtual realdp f(const Variable &x) const;

		virtual Vector &EucGrad(const Variable &x, Vector *result) const;
		virtual Vector &EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const;

		Vector B;
		Vector D;
		integer n;
		integer p;
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of STIEBROCKETT_H */
