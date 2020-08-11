/*
This file defines the class for the Grassmann manifold \Gr(p, n) = \{[X]| X^T X = I_p, X \in R^{n \times p}\}
and [X] = \{XO | O^T O = I_p, O \in R^{p \times p}\}.
It defines the common properties and features of the manifold.

Manifold --> Grassmann

---- WH
*/

#ifndef GRASSMANN_H
#define GRASSMANN_H

#include "Manifolds/Manifold.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class Grassmann : public Manifold{
	public:
		/*Construct the Grassmann manifold: Gr(p, n) and set up default parameters*/
		Grassmann(integer n, integer p);

		/*Delete EMPTYINTR and EMPTYEXTR*/
		virtual ~Grassmann(void);
        
        /*Randomly generate a point on the Grassmann manifold.*/
        virtual Variable RandominManifold(void) const;

		/*Compute the qf retraction defined in [AMS2008, (4.8)].
			[AMS2008]P.-A. Absil, R. Mahony, and R. Sepulchre. Optimization algorithms on matrix manifolds.
			Princeton University Press, Princeton, NJ, 2008.*/
		virtual Variable &Retraction(const Variable &x, const Vector &etax, Variable *result) const;

		/*This cotangent vector is used in the RBFGS defined in [RW2012]
			[RW2012]: W. Ring and B. Wirth. Optimization methods on Riemannian manifolds and their application to shape space.
			SIAM Journal on Optimization, 22(2):596?27, January 2012. */
		virtual Vector &coTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const;

		/*The vector transport by differentiated the retraction of the qf retraction */
		virtual Vector &DiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/*Obtain beta = \|etax\| / \|\mathcal{T}_{R_etax} etax\|
		beta has computed in "DiffRetraction". It is not necessary to recompute it in this function. */
		virtual realdp Beta(const Variable &x, const Vector &etax) const;

		/*Call a member function "ObtainIntrHHR" or "ObtainIntrSquare" based on member variable "retraction". */
		virtual Vector &ObtainIntr(const Variable &x, const Vector &etax, Vector *result) const;

		/*Call a member function "ObtainExtrHHR" or "ObtainExtrSquare" based on member variable "retraction". */
		virtual Vector &ObtainExtr(const Variable &x, const Vector &intretax, Vector *result) const;

		/*ExtrProjection is given by: result = v - x x^T v.*/
		virtual Vector &ExtrProjection(const Variable &x, const Vector &etax, Vector *result) const;

		/*Check whether all the parameters are legal or not.*/
		virtual void CheckParams(void) const;

		/*the Riemannian gradient is obtained by projecting the Euclidean onto the tangent space of x.*/
		virtual Vector &EucGradToGrad(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const;

		/*From the Euclidean action of the Hessian exix:=D grad f[etax] to the Riemannian action of the Hessian:
		xix:=\nabla_etax Rgrad f(x) = (I - x x^T) (exix - etax x^T Egradf)*/
		virtual Vector &EucHvToHv(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const;

	protected:
        mutable integer n; /*The number of row*/
		mutable integer p; /*The number of column*/
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of GRASSMANN_H */
