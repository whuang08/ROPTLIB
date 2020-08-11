/*
This file defines the class for the complex Stiefel manifold \CSt(p, n) = \{X \in C^{n \times p} | X^H X = I_p\}
It defines the common properties and features of the manifold.

Manifold --> CStiefel

---- WH
*/

#ifndef CSTIEFEL_H
#define CSTIEFEL_H

#include "Manifolds/Manifold.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	/*Note that not all metrics, retractions and vector transports have been done.*/

	/* Riemannian Metric for the Stiefel manifold:
	Eucldean: g_x(etax, xix) = Re(\trace(etax^H xix));
	Canonical: g_x(etax, xix) = Re(\trace(etax^H (I_n - x x^H / 2) xix)); */
	enum CStieMetric{ CSTIE_EUCLIDEAN, CSTIE_CANONICAL, CSTIEMETRICLENGTH };

	/*Retraction for the Stiefel manifold
	Complex versions of the retractions in Stiefel.h */
	enum CStieRetraction{ CSTIE_QF, CSTIE_POLAR, CSTIE_EXP, CSTIERETRACTIONLENGTH };

	/*Vector transport for the Stiefel manifold
	Complex versions of the vector transports in Stiefel.h */
	enum CStieVectorTransport{ CSTIE_PARALLELIZATION, CSTIE_PARALLELTRANSLATION, CSTIE_PROJECTION, CSTIEVECTORTRANSPORTLENGTH };

	class CStiefel : public Manifold{
	public:
		/*Construct the Stiefel manifold: CSt(p, n) and set up default parameters*/
		CStiefel(integer n, integer p);

		/*Delete EMPTYINTR and EMPTYEXTR*/
		virtual ~CStiefel(void);

		/* choose Euclidean metric, qf, parallelization and intrinsic representation*/
		virtual void ChooseParamsSet1(void);

		/* choose Euclidean metric,  qf retraction, vector transport by projection and extrinsic representation
		TODO */
		virtual void ChooseParamsSet2(void);

		/* choose Euclidean metric,  polar retraction, vector transport by parallelization and intrinsic representation
		TODO */
		virtual void ChooseParamsSet3(void);

        /* choose Euclidean metric,  polar retraction, vector transport by projection and extrinsic representation
        TODO */
        virtual void ChooseParamsSet4(void);

		/*Euclidean metric*/
        virtual realdp Metric(const Variable &x, const Vector &etax, const Vector &xix) const;

        /*Randomly generate a point in the manifold.*/
        virtual Variable RandominManifold(void) const;

		/*Call a member function "IntrProjection" or "ExtrProjection" based on member variable "IsIntrApproach"*/
        virtual Vector &Projection(const Variable &x, const Vector &etax, Vector *result) const;

		/*Call a member function "qfRetraction" or "ConRetraction" based on member variable "retraction". */
        virtual Variable &Retraction(const Variable &x, const Vector &etax, Variable *result) const;

        /*Compute the inverse of the retraction: result = R_x^{-1}(y).*/
        virtual Vector &InvRetraction(const Variable &x, const Variable &y, Vector *result) const;

		/*Call a member function "qfcoTangentVector" or "ConcoTangentVector" based on member variable "retraction". */
        virtual Vector &coTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const;

		/*Call a member function "DiffqfRetraction" or "DiffConRetraction" based on member variable "retraction". */
        virtual Vector &DiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/*Obtain beta = \|etax\| / \|\mathcal{T}_{R_etax} etax\|
		beta has computed in "DiffRetraction". It is not necessary to recompute it in this function. */
        virtual realdp Beta(const Variable &x, const Vector &etax) const;

		/*The implementation of vector transport by parallization is identity under the intrinsic representation.
		If one needs to use the idea in [HGA2015, Section 4.1], then Manifold::LCVectorTransport is called instead.
		[HGA2015]:Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
		SIAM Journal on Optimization, 25(3):1660?685,2015. */
        virtual Vector &VectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const;

		/*The implementation of inverse vector transport by parallization is identity under the intrinsic representation.
		If one needs to use the idea in [HGA2015, Section 4.1], then Manifold::LCInverseVectorTransport is called instead.
		[HGA2015]:Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
		SIAM Journal on Optimization, 25(3):1660?685,2015. */
        virtual Vector &InverseVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const;

		/*The implementation of inverse vector transport by parallization is identity under the intrinsic representation.
		Therefore, Manifold::HInvTran can be used directly. If one needs to use the idea in [HGA2015, Section 4.1], then
		Manifold::LCHInvTran is called instead.
		[HGA2015]:Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
		SIAM Journal on Optimization, 25(3):1660?685,2015. */
        virtual LinearOPE &HInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, integer start, integer end, LinearOPE *result) const;

		/*The implementation of vector transport by parallization is identity under the intrinsic representation.
		Therefore, Manifold::TranH can be used directly. If one needs to use the idea in [HGA2015, Section 4.1],
		then Manifold::LCTranH is called instead.
		[HGA2015]:Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
		SIAM Journal on Optimization, 25(3):1660?685,2015. */
        virtual LinearOPE &TranH(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, integer start, integer end, LinearOPE *result) const;

		/*The implementation of vector transport by parallization and inverse vector transport by parallelization are identity under
		the intrinsic representation. Therefore, Manifold::TranHInvTran can be used directly. If one needs to use the idea in [HGA2015, Section 4.1],
		then Manifold::LCTranHInvTran is called instead.
		[HGA2015]:Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
		SIAM Journal on Optimization, 25(3):1660?685,2015. */
        virtual LinearOPE &TranHInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, LinearOPE *result) const;

		/*Call a member function "ObtainIntrHHR" based on member variable "retraction". */
        virtual Vector &ObtainIntr(const Variable &x, const Vector &etax, Vector *result) const;

		/*Call a member function "ObtainExtrHHR" based on member variable "retraction". */
        virtual Vector &ObtainExtr(const Variable &x, const Vector &intretax, Vector *result) const;

        /*Compute the intrinsic representation of a vector that is orthogonal to the tangent space at x*/
        virtual Vector &ObtainNorVerIntr(const Variable &x, const Vector &etax, Vector *result) const;

        /*Compute the extrinsic representation of a vector that is orthogonal to the tangent space at x */
        virtual Vector &ObtainNorVerExtr(const Variable &x, const Vector &intretax, Vector *result) const;

		/*IntrProjection is identity, i.e., result <-- v*/
        virtual Vector &IntrProjection(const Variable &x, const Vector &etax, Vector *result) const;

		/*ExtrProjection is given in e.g., [Hua2013, (10.2.5)], i.e., result = v - x sym(x^T v), where sym(M) = (M + M^T)/2.
			[Hua2013]:W. Huang. Optimization algorithms on Riemannian manifolds with applications.
			PhD thesis, Florida State University, Department of Mathematics, 2013.*/
        virtual Vector &ExtrProjection(const Variable &x, const Vector &etax, Vector *result) const;

		/*Check whether all the parameters are legal or not.*/
		virtual void CheckParams(void) const;

		/* Extrinsic representation is used. When the metric is Euclidean, the Riemannian gradient is obtained by projecting the Euclidean
		onto the tangent space of x.*/
        virtual Vector &EucGradToGrad(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const;

		/*When the metric is Euclidean, the Riemannian action of the Hessian is obtained by
		Hess f(x)[etax] = P_x(D grad f(x) [etax]).*/
        virtual Vector &EucHvToHv(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const;

		/*PARAMSMAP is defined in "def.h" and it is a map from string to realdp, i.e., std::map<std::string, realdp> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);

	protected:
		mutable integer n; /*The number of row*/
		mutable integer p; /*The number of column*/
		CStieMetric metric; /*Riemannian metric*/
		CStieRetraction retraction; /*The used retraction*/
		CStieVectorTransport VecTran; /*The used vector transport*/

		/*Householder transformations are used to obtain the intrinsic representation of etax*/
		virtual Vector &ObtainIntrHHR(const Variable &x, const Vector &etax, Vector *result) const;

		/*Householder transformations are used to obtain the extrinsic representation of intretax*/
		virtual Vector &ObtainExtrHHR(const Variable &x, const Vector &intretax, Vector *result) const;

		/*qf retraction defined in [AMS2008, (4.8)]
			[AMS2008]P.-A. Absil, R. Mahony, and R. Sepulchre. Optimization algorithms on matrix manifolds.
			Princeton University Press, Princeton, NJ, 2008.*/
		virtual Variable &qfRetraction(const Variable &x, const Vector &etax, Variable *result) const;

		/*the cotangent vector for the qf retraction in [Hua2013, Section 10.2.4]
			[Hua2013]:W. Huang. Optimization algorithms on Riemannian manifolds with applications.
			PhD thesis, Florida State University, Department of Mathematics, 2013.*/
		virtual Vector &qfcoTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const;

		/*the vector transport by differentiated retraction, see [AMS2008, Example 8.1.5]
			[AMS2008]P.-A. Absil, R. Mahony, and R. Sepulchre. Optimization algorithms on matrix manifolds.
			Princeton University Press, Princeton, NJ, 2008.*/
		virtual Vector &DiffqfRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/*Polar retraction defined in [AMS2008, (4.7)]
		[AMS2008]P.-A. Absil, R. Mahony, and R. Sepulchre. Optimization algorithms on matrix manifolds.
		Princeton University Press, Princeton, NJ, 2008.*/
		virtual Variable &PolarRetraction(const Variable &x, const Vector &etax, Variable *result) const;

        /*Compute the inverse of the retraction: result = R_x^{-1}(y).*/
        virtual Vector &InvPolarRetraction(const Variable &x, const Variable &y, Vector *result) const;

		/*the cotangent vector for the polar retraction in [Hua2013, Section 10.2.4]
		[Hua2013]:W. Huang. Optimization algorithms on Riemannian manifolds with applications.
		PhD thesis, Florida State University, Department of Mathematics, 2013.*/
		virtual Vector &PolarcoTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const;

		/*the vector transport by differentiated retraction, see [HGA2015, Section 4]
		[HGA2015]:W. Huang, K. A. Gallivan, P.-A. Absil, A BROYDEN CLASS OF QUASI-NEWTON METHODS FOR RIEMANNIAN OPTIMIZATION
		SIAM Journal on Optimization, 25(3), 1660-1685, 2015.*/
		virtual Vector &DiffPolarRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir = false) const;
	};
}; /*end of ROPTLIB namespace*/

#endif
