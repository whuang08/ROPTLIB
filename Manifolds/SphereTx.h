/*
This file defines the class for the unit sphere in tangent space of a given manifold at a point of the manifold.

Manifold --> SphereTx

---- WH
*/

#ifndef SPHERETX_H
#define SPHERETX_H

#include "Manifolds/Manifold.h"
#include "Manifolds/Element.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class SphereTx : public Manifold{
	public:
		/*Construct the sphere for the unit sphere in tangent space of a given manifold  at a point of the manifold. */
		SphereTx(Manifold *inmani, Variable *inroot);

		virtual ~SphereTx();

		/*Randomly generate on the unit sphere in the tangent space of the given manifold at a point of the manifold.*/
		Variable RandominManifold() const;

		/*The Riemannian metric is the same as the metric of the input manifold.*/
		virtual realdp Metric(const Variable &x, const Vector &etax, const Vector &xix) const;

		/*etax is in the ambient space. This function projects etax onto the tangent space of x, i.e., result = P_{T_x M} v*/
		virtual Vector &Projection(const Variable &x, const Vector &etax, Vector *result) const;
        
        /*etax is in the ambient space. This function projects etax onto the tangent space of x, i.e., result = P_{T_x M} v*/
        virtual Vector &ExtrProjection(const Variable &x, const Vector &etax, Vector *result) const;
        
        /*Compute result = scalar * etax;*/
        virtual Vector &ScalarTimesVector(const Variable &x, const realdp &scalar, const Vector &etax, Vector *result) const;
        
        /*Compute result = scalar * etax + xix;*/
        virtual Vector &ScalarVectorAddVector(const Variable &x, const realdp &scalar, const Vector &etax, const Vector &xix, Vector *result) const;
        
        /*Compute result = scalar1 * etax + scalar2 * xix;*/
        virtual Vector &VectorLinearCombination(const Variable &x, realdp scalar1, const Vector &etax, realdp scalar2, const Vector &xix, Vector *result) const;

		/*Exponential mapping is used*/
        virtual Variable &Retraction(const Variable &x, const Vector &etax, Variable *result) const;

		/*This is not done yet*/
        virtual Vector &coTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const;

		/*This is not done yet*/
        virtual Vector &DiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/*Return 1; This manifold uses exponential mapping and parallel translation. Therefore, using beta = 1 satisfies the locking condition.*/
        virtual realdp Beta(const Variable &x, const Vector &etax) const;

		/*Parallel translation*/
        virtual Vector &VectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const;

		/*Inverse parallel translation*/
        virtual Vector &InverseVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const;

		/*Check whether all the parameters are legal or not.*/
		virtual void CheckParams(void) const;

		/*When the metric is Euclidean, the Riemannian gradient is obtained by projecting the Euclidean
		onto the tangent space of x.*/
        virtual Vector &EucGradToGrad(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const;

		/*The Riemannian action of the Hessian is obtained by Hess f(x)[etax] = P_x(D grad f(x) [etax]).*/
        virtual Vector &EucHvToHv(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const;
	private:
		Manifold *mani;
		Variable root;
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of SPHERETX_H */
