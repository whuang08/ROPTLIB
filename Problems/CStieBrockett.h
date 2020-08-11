/*
This file defines the class for the problem
min_X tr(X^H B X D), where B is a Hermitian matrix, D is a diagonal real matrix and X \in CSt(p, n).

Problem --> CStieBrockett

---- WH
*/

#ifndef CSTIEBROCKETT_H
#define CSTIEBROCKETT_H

#include "Manifolds/CStiefel.h"
#include "Problems/Problem.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class CStieBrockett : public Problem{
	public:
		CStieBrockett(Vector inB, Vector inD);
		virtual ~CStieBrockett();
		virtual realdp f(const Variable &x) const;

		virtual Vector &EucGrad(const Variable &x, Vector *result) const;
		virtual Vector &EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const;

		Vector B;
		Vector D;
		integer n;
		integer p;
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of CSTIEBROCKETT_H */
