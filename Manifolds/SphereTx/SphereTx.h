/*
This file defines the class for the unit sphere in tangent space of a given manifold at a point of the manifold.

Manifold --> SphereTx

---- WH
*/

#ifndef SPHERETX_H
#define SPHERETX_H

#include "Manifolds/Manifold.h"
//#include "Manifolds/SphereTx/SphereTxVariable.h"
//#include "Manifolds/SphereTx/SphereTxVector.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class SphereTx : public Manifold{
	public:
		/*Construct the sphere for the unit sphere in tangent space of a given manifold  at a point of the manifold. */
		SphereTx(Manifold *inmani, Variable *inroot);

		/*Delete EMPTYINTR and EMPTYEXTR*/
		virtual ~SphereTx();

		/*Randomly generate on the unit sphere in the tangent space of the given manifold at a point of the manifold.
		The returned variable need be deleted by users. Otherwise, it will have memory leakage.*/
		Variable *RandominManifold();

		/*The Riemannian metric is the same as the metric of the input manifold.*/
		virtual double Metric(Variable *x, Vector *etax, Vector *xix) const;

		/*etax is in the ambient space. This function projects etax onto the tangent space of x, i.e., result = P_{T_x M} v*/
		virtual void Projection(Variable *x, Vector *v, Vector *result) const;

		/*Exponential mapping is used*/
		virtual void Retraction(Variable *x, Vector *etax, Variable *result, double stepsize) const;

		/*This is not done yet*/
		virtual void coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*This is not done yet*/
		virtual void DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/*Return 1; This manifold uses exponential mapping and parallel translation. Therefore, using beta = 1 satisfies the locking condition.*/
		virtual double Beta(Variable *x, Vector *etax) const;

		/*Parallel translation*/
		virtual void VectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const;

		/*Inverse parallel translation*/
		virtual void InverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*This is not done yet*/
		virtual void ObtainIntr(Variable *x, Vector *etax, Vector *result) const;

		/*This is not done yet*/
		virtual void ObtainExtr(Variable *x, Vector *intretax, Vector *result) const;

		/*result <-- v*/
		virtual void IntrProjection(Variable *x, Vector *v, Vector *result) const;

		/*etax is in the ambient space. This function projects etax onto the tangent space of x, i.e., result = P_{T_x M} v*/
		virtual void ExtrProjection(Variable *x, Vector *v, Vector *result) const;

		/*Check whether all the parameters are legal or not.*/
		virtual void CheckParams(void) const;

		/*When the metric is Euclidean, the Riemannian gradient is obtained by projecting the Euclidean
		onto the tangent space of x.*/
		virtual void EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const;

		/*The Riemannian action of the Hessian is obtained by Hess f(x)[etax] = P_x(D grad f(x) [etax]).*/
		virtual void EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const;
	private:
		Manifold *mani;
		Variable *root;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of L2SPHERE_H
