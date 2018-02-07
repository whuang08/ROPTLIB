/*
This file defines the blind deconvolution problem
min_{X \in R_r^{m by n}} diag(F B X (F C)^T), 
where R_r{m by n} is the set of m by n matrices with rank r.

Problem --> LRBlindDeconvolution

---- WH
*/

#ifndef LRBLINDDECONVOLUTION_H
#define LRBLINDDECONVOLUTION_H

#include "Manifolds/Stiefel/Stiefel.h"
#include "Manifolds/Stiefel/StieVariable.h"
#include "Manifolds/Stiefel/StieVector.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Others/def.h"
#include "Manifolds/Element.h"
#include "Manifolds/ProductElement.h"
#include "Manifolds/ProductManifold.h"
#include "Manifolds/LowRank/LowRank.h"
#include "Others/fftw/fftw3.h"
#include "Others/SparseBLAS/blas_sparse.h"

#ifdef ROPTLIB_WITH_FFTW

/*Define the namespace*/
namespace ROPTLIB{

	class LRBlindDeconvolution : public Problem{
	public:
		/*y \in C^L, B \in C^{L \times K}, C \in C^{L \times N}, */
		LRBlindDeconvolution(double *iny, double *inB, integer innzmaxB, int *inirB, int *injcB, bool inisBsparse, 
			double *inC, integer innzmaxC, int *inirC, int *injcC, bool inisCsparse, integer inL, integer inK, integer inN, integer inr);
		virtual ~LRBlindDeconvolution();
		virtual double f(Variable *x) const;

		virtual void RieGrad(Variable *x, Vector *gf) const;
		virtual void RieHessianEta(Variable *x, Vector *etax, Vector *xix) const;
		void fftwrapper(integer n, integer r, fftw_complex *in, fftw_complex *out, int sign) const;
		void BtimesU(double *U, double *result) const;
		void CtimesV(double *V, double *result) const;
		void BTtimesM(double *M, double *result) const;
		void CTtimesM(double *V, double *result) const;

		double *y;
		
		double *B;
		integer nzmaxB;
		int *irB;
		int *jcB;
		bool isBsparse;
		blas_sparse_matrix sB;

		double *C;
		integer nzmaxC;
		int *irC;
		int *jcC;
		bool isCsparse;
		blas_sparse_matrix sC;

		integer L;
		integer K;
		integer N;
		integer r;
		mutable fftw_plan p;
		unsigned int flags;
	};
}; /*end of ROPTLIB namespace*/
#endif
#endif
