/*
This file defines the Phase Retrieval problem
min_{U \in C_*^{n \times p} / U_p} \|b - diag(Z U U^* Z^*)\|_2^2 + \kappa \trace(X),
where n = n1 * n2, b \in R^m, m = n1 * n2 * l

Problem --> CSFRQPhaseRetrieval

---- WH
*/

/*
 TODO, the efficiency can be improved.
 *), try not to use submatrix
 *), set fft_plan only once in FFT2D
*/

#ifndef CSFRQPHASERETRIEVAL_H
#define CSFRQPHASERETRIEVAL_H

#include "Problems/Problem.h"
#include "Others/def.h"
#include "Manifolds/Element.h"
#include "Manifolds/CSymFixedRankQ.h"

#ifdef ROPTLIB_WITH_FFTW

/*Define the namespace*/
namespace ROPTLIB{

	class CSFRQPhaseRetrieval : public Problem{
	public:
		/*b \in R^m, m = n1*n2*l, masks \in C^{n \times l}*/
		CSFRQPhaseRetrieval(Vector inb, Vector inmasks, realdp inkappa, integer inn1, integer inn2, integer inl, integer inr);

		virtual ~CSFRQPhaseRetrieval();
		virtual realdp f(const Variable &x) const;

		virtual Vector &EucGrad(const Variable &x, Vector *result) const;
        
        virtual Vector &EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const;
        
        Vector b;
        Vector masks;
        
		realdp kappa;

		mutable integer n1;
        mutable integer n2;
        mutable integer n;
		mutable integer l;
		mutable integer r;
		mutable integer m;
	};
}; /*end of ROPTLIB namespace*/
#endif
#endif /*CSFRQPHASERETRIEVAL_H*/
