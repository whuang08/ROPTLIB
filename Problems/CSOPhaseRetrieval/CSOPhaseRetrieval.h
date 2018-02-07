/*
This file defines the Phase Retrieval problem
min_{U \in C_*^{n \times p} / O_p} \|b - diag(Z U U^* Z^*)\|_2^2 + \kappa \trace(X), 

Problem --> CSOPhaseRetrieval

---- WH
*/

#ifndef CSOPHASERETRIEVAL_H
#define CSOPHASERETRIEVAL_H

#include "Manifolds/CpxNStQOrth/CpxNStQOrth.h"
#include "Manifolds/CpxNStQOrth/CSOVariable.h"
#include "Manifolds/CpxNStQOrth/CSOVector.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Others/def.h"
#include "Others/fftw/fftw3.h"

#ifdef ROPTLIB_WITH_FFTW

/*Define the namespace*/
namespace ROPTLIB{

	class CSOPhaseRetrieval : public Problem{
	public:
		/*b \in R^m, m = n1*n2*l, masks \in C^{n1 \times n2 \times l}*/
		CSOPhaseRetrieval(double *inb, double *inmasks, double inkappa, integer inn1, integer inn2, integer inl, integer inr);

		virtual ~CSOPhaseRetrieval();
		virtual double f(Variable *x) const;

		virtual void EucGrad(Variable *x, Vector *gf) const;

		void fftwrapper(integer inn1, integer inn2, fftw_complex *in, fftw_complex *out, int sign) const;

		double *b;
		double *masks;
		double kappa;

		integer n1;
		integer n2;
		integer l;
		integer r;
		integer m;
		mutable fftw_plan p;
		unsigned int flags;
	};
}; /*end of ROPTLIB namespace*/
#endif
#endif
