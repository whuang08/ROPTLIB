/*
This file defines the class for the problem: 
Computing the karcher mean of SPD manifold, i.e.,
min_{X \in SPD} (1 / 2 / num) sum_{i=1}^{num} dist^2(X, A_i)
where A_i \in SPD and dist(A, B) = \|logm(A^{-1/2} B A^{-1/2})\|_F.

Problem --> SPDKarcherMean

---- WH
*/

#ifndef SPDKARCHERMEAN_H
#define SPDKARCHERMEAN_H

#include "Manifolds/SPDManifold.h"
#include "Problems/Problem.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class SPDKarcherMean : public Problem{
	public:
		/*The sample SPD matrices Ai are stored as their Cholesky matrix, i.e., Ai = Li Li^T.
		inLs is an array of Vectors with length innum.*/
		SPDKarcherMean(Vector inLs, integer inn, integer innum);
		virtual ~SPDKarcherMean();
		virtual realdp f(const Variable &x) const;

		virtual Vector &EucGrad(const Variable &x, Vector *result) const;

		virtual Vector &EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const;

		Vector Ls;
		integer n;
		integer num;
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of STIEBROCKETT_H */
