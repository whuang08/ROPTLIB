/*
This file defines the class for the the complex quotient manifold C_*^{m \times r} \times C_*^{n \times r} / GL(r), 
which is equivalent to the complex fixed-rank manifold C_r^{m \times n},
where C_*^{s \times t} denotes the complex noncompact Stiefel manifold and GL(r) denotes the complex general group. 

Riemannian metric is different from the Euclidean metric, it is \trace(dG1^* dG2 (H^* H) + dH1^* dH2 (G^* G))

Manifold --> ProductManifold --> CFixedRank2Factors

---- WH
*/

#ifndef CFIXEDRANK2FACTORS_H
#define CFIXEDRANK2FACTORS_H

#include "Manifolds/ProductManifold.h"
#include "Manifolds/Euclidean/Euclidean.h"
#include "Manifolds/CFixedRank2Factors/CFR2Variable.h"
#include "Manifolds/CFixedRank2Factors/CFR2Vector.h"
#include "Manifolds/CpxNStQOrth/CpxNStQOrth.h"
#include "Others/MyMatrix.h"

/*Define the namespace*/
namespace ROPTLIB{

	class CFixedRank2Factors : public ProductManifold{
	public:
		/*Construct the complex quotient manifold C_*^{m \times r} \times C_*^{n \times r} / GL(r),
		which is equivalent to the complex fixed-rank manifold C_r^{m \times n}.*/
		CFixedRank2Factors(integer m, integer n, integer r);

		/*Delete the manifold.*/
		~CFixedRank2Factors(void);

		/*Intrinsic representation has size (m + n - r) r*/
		virtual void ObtainIntr(Variable *x, Vector *etax, Vector *result) const;

		/*Compute the extrinsic approach given the intrinsic representation*/
		virtual void ObtainExtr(Variable *x, Vector *intretax, Vector *result) const;

		/*Perform the default retraction of each manifold component*/
		virtual void Retraction(Variable *x, Vector *etax, Variable *result, double stepsize) const;

		/*Compute the tangent vector result satisfying
		g_y(\mathcal{T}_{R_etax}(xix), xiy) = g_x(xix, result) for all xix \in T_x M,
		where y = R_x(etax), xiy \in T_y M. (TODO)*/
		virtual void coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*Perform the vector transport by differentiated retraction.*/
		virtual void DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const;
		
		/*etax is in the total space C_*^{m \times r} \times C_*^{n \times r}. This function projects etax onto the 
		tangent space of x, i.e., result = P_{T_x M} v, where P is based on the selected Riemannian metric*/
		virtual void ExtrProjection(Variable *x, Vector *v, Vector *result) const;

		/*the Riemannian gradient is obtained by projecting the Euclidean onto the tangent space of x.*/
		virtual void EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const;

		/*The Riemannian action of the Hessian is obtained by Hess f(x)[etax] = P_x(D grad f(x) [etax]).*/
		virtual void EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const;

	protected:
		mutable integer m; /*the number of row*/
		mutable integer n; /*the number of column*/
		mutable integer r; /*the rank of the matrix*/
	};
}; /*end of ROPTLIB namespace*/
#endif // end of LOWRANK_H
