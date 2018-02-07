/*
This file defines the blind deconvolution problem
min_{h \in R^K, m \in R^N}} \|y - diag(F B h m^T (F C)^T)\|_2^2.

Problem --> EucBlindDeconvolution

---- WH
*/

#ifndef EUCBLINDDECONVOLUTION_H
#define EUCBLINDDECONVOLUTION_H

#include "Manifolds/Euclidean/Euclidean.h"
#include "Manifolds/Euclidean/EucVariable.h"
#include "Manifolds/Euclidean/EucVector.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Others/def.h"
#include "Manifolds/Element.h"
#include "Others/fftw/fftw3.h"
#include "Others/SparseBLAS/blas_sparse.h"

#ifdef ROPTLIB_WITH_FFTW

/*Define the namespace*/
namespace ROPTLIB{

	class EucBlindDeconvolution : public Problem{
	public:
		/*y \in C^L, B \in C^{L \times K}, C \in C^{L \times N}, */
		EucBlindDeconvolution(double *iny, double *inB, integer innzmaxB, int *inirB, int *injcB, bool inisBsparse,
			double *inC, integer innzmaxC, int *inirC, int *injcC, bool inisCsparse, integer inL, integer inK, integer inN, integer inr);
		virtual ~EucBlindDeconvolution();
		virtual double f(Variable *x) const;

		virtual void EucGrad(Variable *x, Vector *gf) const;
		virtual void EucHessianEta(Variable *x, Vector *etax, Vector *xix) const;
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
