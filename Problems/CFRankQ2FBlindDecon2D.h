/*
This file defines the blind deconvolution problem for 2D n1-by-n2 problems
min_{(h, m) \in C_*^N \times C_*^N / C_*} \|y - diag(F B h (\bar{F} H C m)^*)\|_2^2 + rho sum_{i=1}^L max( (L |b_i^* h|^2 \|m\|_2^2 ) / (8 d^2 mu^2) - 1, 0)^2,
where y \in C^L, B \in C^{L \times N}, C \in C^{L \times N}, H \in C^{L \times L},  L = N.
This is for 2D problems. F is the Kronecker Product of two 1D DFT matrix.
See the optimization problem in
    Wen Huang and Paul Hand, Blind Deconvolution by a Steepest Descent Algorithm on a Quotient Manifold, SIIMS, 2018

Note that when this is used for blind deconvolution in imaging denoising, B is the
is a sparse diagonal matrix where the sparse entries correspond to the support of the kernel;
and H is the Haar wavelet transform matrix and C is a sparse diagonal matrix where the sparse
entries correspond to the columns of H that are chosen.
Suppose the image is n1 by n2 resolution, then L = N = n1 * n2;

Problem --> CFRankQ2FBlindDecon2D

---- WH
*/

#ifndef CFRANKQ2FBLINDDECON2D_H
#define CFRANKQ2FBLINDDECON2D_H

#include "Problems/Problem.h"
#include "Others/def.h"
#include "Manifolds/Element.h"
#include "Manifolds/CFixedRankQ2F.h"
#undef abs
#include "Others/wavelet/wavelet.h"

#ifdef ROPTLIB_WITH_FFTW

/*Define the namespace*/
namespace ROPTLIB{

	class CFRankQ2FBlindDecon2D : public Problem{
	public:
		/*y \in C^L, B \in C^{L \times N}, C \in C^{L \times N},
		if C == nullptr, then haar wavelet transform is used.*/
		CFRankQ2FBlindDecon2D(Vector iny, SparseMatrix &inB, SparseMatrix &inC, integer inn1, integer inn2, integer inr, realdp inrho, realdp ind, realdp inmu);
		virtual ~CFRankQ2FBlindDecon2D();
    
        virtual realdp f(const Variable &x) const;

        virtual Vector &EucGrad(const Variable &x, Vector *result) const;

        virtual Vector &EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const;
        
		Vector y;
        SparseMatrix *Bptr;
        SparseMatrix *Cptr;

		mutable integer L;
		mutable integer n1;
		mutable integer n2;
		mutable integer r;
        
		realdp rho;
		realdp d;
		realdp mu;
	};
}; /*end of ROPTLIB namespace*/
#endif

#endif
