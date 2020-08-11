/*
This file defines the class for the problem
min_{X \in R_r^{m by n}} \|X - A\|_W^2, 
where R_r{m by n} is the set of m by n matrices with rank r, A is a given m by n matrix, W is a mn by mn symmetric
positive definite matrix and \|M\|_W^2 is the W weighted norm, i.e., \|M\|_W^2 = vec(M)^T W vec(M) and vec(M) is the vector given by
vectorizing the matrix M.
This problem is used in [ZHGVA2015]
	[ZHGVA2015]: Guifang Zhou, Wen Huang, Kyle A. Gallivan, Paul Van Dooren, Pierre-Antoine Absil. Rank-Constrained Optimization: A Riemannian Manifold Approach,
		In Proceeding of European Symposium on Artificial Neural Networks, Computational Intelligence and Machine Learning, 2015.

Problem --> FRankWeightApprox

---- WH
*/

#ifndef FRANKEWEIGHTAPPROX_H
#define FRANKEWEIGHTAPPROX_H

#include "Problems/Problem.h"
#include "Others/def.h"
#include "Manifolds/FixedRankE.h"

/*Define the namespace*/
namespace ROPTLIB{

	class FRankEWeightApprox : public Problem{
	public:
		FRankEWeightApprox(Vector inA, Vector inW, integer inm, integer inn, integer inr);
		virtual ~FRankEWeightApprox();
		virtual realdp f(const Variable &x) const;

		virtual Vector &EucGrad(const Variable &x, Vector *result) const;
		virtual Vector &EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const;

		Vector A;
		Vector W;
		integer m;
		integer n;
		integer r;
	};
}; /*end of ROPTLIB namespace*/
#endif
