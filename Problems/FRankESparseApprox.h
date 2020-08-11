/*
This file defines the class for the problem
min_{X \in R_r^{m by n}} \|X - A\|_F^2 + lambda \|X\|_1,
where R_r{m by n} is the set of m by n matrices with rank r, A is a given m by n matrix.

Problem --> FRankESparseApprox

---- WH
*/

#ifndef FRANKESPARSEAPPROX_H
#define FRANKESPARSEAPPROX_H

#include "Problems/Problem.h"
#include "Others/def.h"
#include "Manifolds/FixedRankE.h"

/*Define the namespace*/
namespace ROPTLIB{

	class FRankESparseApprox : public Problem{
	public:
		FRankESparseApprox(Vector inA, realdp inlambda, integer inm, integer inn, integer inr, integer inlengthW);
		virtual ~FRankESparseApprox();
		virtual realdp f(const Variable &x) const;

		virtual Vector &EucGrad(const Variable &x, Vector *result) const;
		virtual Vector &EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const;
        
        virtual Vector &ProxW(const Vector &x, const Vector &Weight, Vector *result) const;
        virtual Vector &CalJW(const Vector &x, const Vector &eta, const Vector &Weight, Vector *result) const;
        
        /*The create the weight matrix*/
        virtual Vector &PreConditioner(const Variable &x, const Vector &eta, Vector *result) const;

		Vector A;
        realdp lambda;
		integer m;
		integer n;
		integer r;
        integer lengthW;
        realdp L;
	};
}; /*end of ROPTLIB namespace*/
#endif
