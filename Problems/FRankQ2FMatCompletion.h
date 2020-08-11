/*
This file defines the class for the problem
min_{X \in R_r^{m by n}} 0.5 \|P_{\Omega}(X) - P_{\Omega}(A)\|_F^2, 
where R_r{m by n} is the set of m by n matrices with rank r, \Omega is a index set, A_{ij}, (i, j) \in \Omega are given.
The fixed rank manifold R_r^{m times n} is represented by a 2 factor quotient manifold

Problem --> FRankQ2FMatCompletion

---- WH
*/

#ifndef FRANKQ2MATCOMPLETION_H
#define FRANKQ2MATCOMPLETION_H

#include "Problems/Problem.h"
#include "Others/def.h"
#include "Manifolds/Element.h"
#include "Manifolds/FixedRankQ2F.h"

/*Define the namespace*/
namespace ROPTLIB{

	class FRankQ2FMatCompletion : public Problem{
	public:
		FRankQ2FMatCompletion(integer *inir, integer *injc, realdp *invals, integer innz, integer inm, integer inn, integer inr);

		virtual ~FRankQ2FMatCompletion();

		/*0.5 \|P_omaga(X) - P_omega(A)\|_F^2*/
		virtual realdp f(const Variable &x) const;

		/*P_omaga(X) - P_omega(A)*/
		virtual Vector &EucGrad(const Variable &x, Vector *result) const;

		/*P_omaga(etax)*/
		virtual Vector &EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const;

        /*compute result = P_Omega(G H^T), where the indices of Omega are stored in inir and injc*/
        void ProjecOmegaGHT(const realdp *G, const realdp *H, integer inm, integer inn, integer inr, integer *inir, integer *injc, integer nz, realdp *result) const;
        
        integer *ir;
        integer *jc;
        realdp *vals;
        mutable integer nz;
		mutable integer m;
		mutable integer n;
		mutable integer r;
	};
}; /*end of ROPTLIB namespace*/
#endif
