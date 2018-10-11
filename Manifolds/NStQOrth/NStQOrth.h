/*
This file defines the class for the manifold R_*^{n \times p} / O_p, where R_*^{n \times p} is a n by p full column rank
matrix and O_p is a p-by-p orthogonal group. This manifold is equivalent to the manifold of the set of symmetric positive 
semidefinite matricies with rank fixed (rank p). The Riemannian metric is
   g_Y(\eta_Y, \xi_Y) = 2 tr(Y^T \eta_Y Y^T \xi_Y + Y^T Y \eta_Y^T \xi_Y)
This Riemannian metric is equivalent to the Euclidean metric on the set of symmetric positive semidefinite matricies with 
rank fixed (rank p)

Manifold --> NStQOrth

---- WH
*/

#ifndef NSTQORTH_H
#define NSTQORTH_H

#include "Manifolds/NStQOrth/NSOVariable.h"
#include "Manifolds/NStQOrth/NSOVector.h"
#include "Manifolds/Manifold.h"
#include "Others/MyMatrix.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	/*Note that not all metrics, retractions and vector transports have been done.*/

	/* Riemannian Metric for the R_*^{n \times p} / O_p manifold:
	EUC: g_x(etax, xix) = 2 \trace(x^T etax x^T xix + x^T x etax^T xix); This metric is equivalent to the Euclidean metric on the symmetric positive semidefinite matrix with rank p
	QEUC, g_x(etax, xix) = \trace(etax^T xix);
	HGZ, g_x(etax, xix) = \trace(x^T x etax^T xix)
	*/
	enum NSOMetric { NSOEUC, NSOQEUC, NSOHGZ, NSOMETRICLENGTH };

	/*Vector transport:
	NSOPARALLELIZATION: vector transport by parallelization
	NSOPROJECTION: vector transport by projection
	*/
	enum NSOVectorTransport { NSOPARALLELIZATION, NSOPROJECTION, NSOVECTORTRANSPORTLENGTH };

	class NStQOrth : public Manifold{
	public:
		/*Construct the manifold C_*^{n \times p} / O_p*/
		NStQOrth(integer n, integer p = 1);

		/*Choose Riemannian metric: g_x(etax, xix) = \trace(x^T x etax^T xix), vector transport by parallelization */
		void ChooseNSOParamsSet1(void);

		/*Choose Riemannian metric: g_x(etax, xix) = 2 \trace(x^T etax x^T xix + x^T x etax^T xix), vector transport by parallelization */
		void ChooseNSOParamsSet2(void);

		/*Choose Riemannian metric: g_x(etax, xix) = \trace(etax^T xix), vector transport by parallelization */
		void ChooseNSOParamsSet3(void);

		/*Choose Riemannian metric: g_x(etax, xix) = \trace(x^T x etax^T xix), vector transport by projection */
		void ChooseNSOParamsSet4(void);

		/*Choose Riemannian metric: g_x(etax, xix) = 2 \trace(x^T etax x^T xix + x^T x etax^T xix), vector transport by projection */
		void ChooseNSOParamsSet5(void);

		/*Choose Riemannian metric: g_x(etax, xix) = \trace(etax^T xix), vector transport by projection */
		void ChooseNSOParamsSet6(void);

		/*Delete EMPTYINTR and EMPTYEXTR*/
		virtual ~NStQOrth(void);

		/*Euclidean metric*/
		virtual double Metric(Variable *x, Vector *etax, Vector *xix) const;

		/*Call a member function "IntrProjection" or "ExtrProjection" based on member variable "IsIntrApproach"*/
		virtual void Projection(Variable *x, Vector *v, Vector *result) const;

		/*Only support intrinsic representataion for etax. The retraction is result = x + etax.*/
		virtual void Retraction(Variable *x, Vector *etax, Variable *result, double stepsize) const;

		/*This is not done yet. */
		virtual void coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*The vector transport by differentiated retraction is the same is the vector tranport by projection for this retraction and manifold.*/
		virtual void DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/*The vector transport*/
		virtual void VectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const;

		/*Only use one*/
		virtual double Beta(Variable *x, Vector *etax) const;

		/*Householder transformations are used to obtain the intrinsic representation of etax*/
		virtual void ObtainIntr(Variable *x, Vector *etax, Vector *result) const;

		/*Householder transformations are used to obtain the extrinsic representation of etax*/
		virtual void ObtainExtr(Variable *x, Vector *intretax, Vector *result) const;

		/*IntrProjection is identity, i.e., result <-- v*/
		virtual void IntrProjection(Variable *x, Vector *v, Vector *result) const;

		/*ExtrProjection*/
		virtual void ExtrProjection(Variable *x, Vector *v, Vector *result) const;

		/*ExtrProjection for the QEUC Riemannian metric*/
		virtual void ExtrProjectionQEUC(Variable *x, Vector *v, Vector *result) const;

		/*ExtrProjection for the QEUC Riemannian metric*/
		virtual void ExtrProjectionEUCorHGZ(Variable *x, Vector *v, Vector *result) const;

		/*Check whether all the parameters are legal or not.*/
		virtual void CheckParams(void) const;

		/**/
		virtual void EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const;

		/**/
		virtual void EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const;

		/*Compute the units vectors in Hourseholder transformations, which can be used in the member functions "ObtainIntr" and "ObtainExtr". */
		static void ComputeHHR(Variable *x);
	protected:

		/*for EUC metric*/
		virtual void EucGradToGradEUC(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const;

		/*for EUC metric*/
		virtual void EucHvToHvEUC(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const;

		/*for EUC metric*/
		virtual void ObtainIntrEUC(Variable *x, Vector *etax, Vector *result) const;

		/*for EUC metric*/
		virtual void ObtainExtrEUC(Variable *x, Vector *intretax, Vector *result) const;

		/*for QEUC metric*/
		virtual void EucGradToGradQEUC(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const;

		/*for QEUC metric*/
		virtual void EucHvToHvQEUC(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const;

		/*for QEUC metric*/
		virtual void ObtainIntrQEUC(Variable *x, Vector *etax, Vector *result) const;

		/*for QEUC metric*/
		virtual void ObtainExtrQEUC(Variable *x, Vector *intretax, Vector *result) const;

		/*for HGZ metric*/
		virtual void EucGradToGradHGZ(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const;

		/*for HGZ metric*/
		virtual void EucHvToHvHGZ(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const;

		/*for HGZ metric*/
		virtual void ObtainIntrHGZ(Variable *x, Vector *etax, Vector *result) const;

		/*for HGZ metric*/
		virtual void ObtainExtrHGZ(Variable *x, Vector *intretax, Vector *result) const;
		
		/*The vector transport by projection*/
		virtual void VectorTransportProj(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const;

		integer n; /*The row*/
		integer p; /*The column*/
		NSOMetric metric;
		NSOVectorTransport VecTran;
	};
}; /*end of ROPTLIB namespace*/

#endif // end of NSTQORTH_H
