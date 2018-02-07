/*
This file defines the blind deconvolution problem
min_{(h, m) \in C_*^K \times C_*^N / C_*} \|y - diag(F B h (\bar{F} C m)^*)\|_2^2,
where $m$ and $C$ is different from the two in EucBlindDeconvolution. m_here = \bar{m_there}, C_here = \bar{C_there}

Problem --> CFR2BlindDeconvolution

---- WH
*/

#ifndef CFR2BLINDDECONVOLUTION_H
#define CFR2BLINDDECONVOLUTION_H

#include "Manifolds/CFixedRank2Factors/CFixedRank2Factors.h"
#include "Manifolds/CFixedRank2Factors/CFR2Variable.h"
#include "Manifolds/CFixedRank2Factors/CFR2Vector.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Others/def.h"
#include "Manifolds/Element.h"
#include "Others/fftw/fftw3.h"
#include "Others/SparseBLAS/blas_sparse.h"
#undef abs
#include "Others/wavelet/wavelet.h"

//#define max(a, b) ((a < b)?b:a)

#ifdef ROPTLIB_WITH_FFTW

/*Define the namespace*/
namespace ROPTLIB{

	class CFR2BlindDeconvolution : public Problem{
	public:
		/*y \in C^L, B \in C^{L \times K}, C \in C^{L \times N},
		if C == nullptr, then haar wavelet transform is used.*/
		CFR2BlindDeconvolution(double *iny, double *inB, integer innzmaxB, int *inirB, int *injcB, bool inisBsparse,
			double *inC, integer innzmaxC, int *inirC, int *injcC, bool inisCsparse, integer inL, integer inK, integer inN, integer inr, double inrho, double ind, double inmu);
		virtual ~CFR2BlindDeconvolution();
		virtual double f(Variable *x) const;

		virtual void EucGrad(Variable *x, Vector *gf) const;
		virtual void EucHessianEta(Variable *x, Vector *etax, Vector *xix) const;
		void fftwrapper(integer n, integer r, fftw_complex *in, fftw_complex *out, int sign) const;
		void BtimesU(double *U, double *result) const;
		void CtimesV(double *V, double *result) const;
		void BHtimesM(double *M, double *result) const;
		void CHtimesM(double *V, double *result) const;

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

		mutable std::vector<double> flagR, flagI;
		mutable std::vector<int> lengthR, lengthI;

		integer L;
		integer K;
		integer N;
		integer r;
		//integer log2L;
		mutable fftw_plan p;
		unsigned int flags;

		double rho;
		double d;
		double mu;
	};
}; /*end of ROPTLIB namespace*/
#endif
#endif
