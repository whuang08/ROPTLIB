/*
This file defines the class for a sparse PCA model on the Stiefel manifold
min_X - tr(X^T B^T B X) + \lambda \|X\|_1, where B is a m-by-n matrix, and X \in St(p, n).
Typically, n > m > p.

Problem --> SPCA

---- WH
*/

#ifndef STIESPCA_H
#define STIESPCA_H

#include "Manifolds/Stiefel.h"
#include "Problems/Problem.h"
#include "Solvers/ManPG.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class StieSPCA : public Problem{
	public:
		StieSPCA(Vector inB, realdp inlambda, integer inn, integer inm, integer inp, integer inlengthW);
		virtual ~StieSPCA();
		virtual realdp f(const Variable &x) const;
        
        virtual Vector &EucGrad(const Variable &x, Vector *result) const;
        virtual Vector &EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const;
        
        virtual Vector &ProxW(const Vector &x, const Vector &Weight, Vector *result) const;
        virtual Vector &CalJW(const Vector &x, const Vector &eta, const Vector &Weight, Vector *result) const;
        
        /*The create the weight matrix*/
        virtual Vector &PreConditioner(const Variable &x, const Vector &eta, Vector *result) const;

		Vector B;
        realdp lambda;
		mutable integer n;
        mutable integer m;
		mutable integer p;
        integer lengthW;
        Vector colnormB;
        realdp L;
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of GRASSRQ_H */
