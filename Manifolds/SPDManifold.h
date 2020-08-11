/*
This file defines the class for the manifold of symmetric positive definite matrices (SPDManifold).
The affine invariant metric is used. Intrinsic representation of tangent vectors are used.
See details in: 
	[YHAG15] Xinru Yuan, Wen Huang, P.-A. Absil, K. A. Gallivan. "A Riemannian Limited-memory BFGS 
	Algorithm for Computing the Matrix Geometric Mean".

Manifold --> SPD

---- WH
*/
#ifndef SPDMANIFOLD_H
#define SPDMANIFOLD_H

#include "Manifolds/Manifold.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

    enum SPDMetric { SPDEUCLIDEAN, SPDAFFINEINVARIANCE, SPDMETRICLENGTH };

    enum SPDRetraction { SPDEXP, SPDSECONDORDER, SPDRETRACTIONLENGTH };

    enum SPDVectorTransport { SPDPARATRAN, SPDVTPARA, SPDVECTORTRANSPORTLENGTH };

	class SPDManifold : public Manifold{
	public:
		/*Construct the SPD manifold*/
		SPDManifold(integer inn);

		/*Delete EMPTYINTR and EMPTYEXTR*/
		virtual ~SPDManifold(void);

		/* choose Affine invariance metric, second order retraction, parallelization and intrinsic representation*/
		virtual void ChooseParamsSet1(void);

		/* choose Euclidean metric, second order retraction, parallelization and intrinsic representation*/
		virtual void ChooseParamsSet2(void);
        
        /* choose Affine invariance metric, second order retraction, parallelization and extrinsic representation*/
        virtual void ChooseParamsSet3(void);
        
        /* choose Euclidean metric, second order retraction, parallelization and extrinsic representation*/
        virtual void ChooseParamsSet4(void);
        
        /*Randomly generate a point on the symmetric positive definite manifold.*/
        virtual Variable RandominManifold(void) const;

		/*Check whether all the parameters are legal or not.*/
		virtual void CheckParams(void) const;
        
        /* metrics */
        virtual realdp Metric(const Variable &x, const Vector &etax, const Vector &xix) const;
        
		/*The second order approximation of the exponential mapping, i.e., result = x + etax + 0.5 etax x^{-1} etax */
		virtual Variable &Retraction(const Variable &x, const Vector &etax, Variable *result) const;
        
        /*The vector transport by differentiated the retraction, which is R_x(etax) = x + etax + 0.5 etax x^{-1} etax.
        Therefore, it is xix + 0.5 xix x^{-1} etax + 0.5 etax x^{-1} xix */
        virtual Vector &DiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/*computes beta = \|etax\| / \|\mathcal{T}_{R_etax} etax\|
		Default: beta <-- 1*/
		virtual realdp Beta(const Variable &x, const Vector &etax) const;

		/*compute the distance between two points on the manifold.
		Default: the distance under affine invariant metric: ||log(x1^{-1/2) x2 x1^{-1/2}||_F*/
        virtual realdp Dist(const Variable &x1, const Variable &x2) const;

		/*TODO*/
		virtual Vector &coTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const;

		/*The orthogonal projection onto the tangent space using the affine invariance metric is: result = (etax + etax^T) / 2*/
		virtual Vector &ExtrProjection(const Variable &x, const Vector &etax, Vector *result) const;

		/*Compute the intrinsic representation of a tagnent vector etax, intreta = upper triangle
		of L^{-1} etax L^{-T}, where x = L L^T. */
		virtual Vector &ObtainIntr(const Variable &x, const Vector &etax, Vector *result) const;

		/*Compute the extrinsic representation of a tagnent vector intretax. Inverse operation of the function ObtainIntr */
		virtual Vector &ObtainExtr(const Variable &x, const Vector &intretax, Vector *result) const;

		/*gf = x * egf * x */
		virtual Vector &EucGradToGrad(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const;

		/*Not done yet. Temporarily use: xix <-- exix*/
		virtual Vector &EucHvToHv(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const;

        /*Computes the vector transport, i.e., result = \mathcal{T}_etax (xix)*/
        virtual Vector &VectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const;
        
	protected:
		integer n; /*The size of the space, i.e., the number of row/column */
		SPDMetric metric; /*Riemannian metric*/
        SPDRetraction retr;
        SPDVectorTransport VecTran;
        
        /*The exponential mapping, i.e., result = x^{1/2} exp(x^{-1/2} etax x^{-1/2}) x^{1/2} */
        virtual Variable &ExpRetraction(const Variable &x, const Vector &etax, Variable *result) const;
        
        /*The second order approximation of the exponential mapping, i.e., result = x + etax + 0.5 etax x^{-1} etax */
        virtual Variable &SecOrdRetraction(const Variable &x, const Vector &etax, Variable *result) const;
        
        /*The vector transport by differentiated the retraction, which is R_x(etax) = x^{1/2} exp(x^{-1/2} etax x^{-1/2}) x^{1/2}. */
        virtual Vector &DiffExpRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir = false) const;
        
        /*The vector transport by differentiated the retraction, which is R_x(etax) = x + etax + 0.5 etax x^{-1} etax.
        Therefore, it is xix + 0.5 xix x^{-1} etax + 0.5 etax x^{-1} xix */
        virtual Vector &DiffSecOrdRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir = false) const;
        
        virtual Vector &ParallelTranslation(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const;
        
        virtual Vector &VectorTransportParallelization(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const;

		/*Compute the intrinsic representation of a tagnent vector etax, intreta = upper triangle
		of L^{-1} etax L^{-T}, where x = L L^T. */
		virtual Vector &ObtainIntrAF(const Variable &x, const Vector &etax, Vector *result) const;

		/*Compute the extrinsic representation of a tagnent vector intretax. Inverse operation of the function ObtainIntrAF */
		virtual Vector &ObtainExtrAF(const Variable &x, const Vector &intretax, Vector *result) const;

		/*gf = x * egf * x */
		virtual Vector &EucGradToGradAF(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const;

		/*Not done yet. Temporarily use: xix <-- exix*/
		virtual Vector &EucHvToHvAF(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const;

		/*compute the distance between two points on the manifold.
		Default: the distance under affine invariant metric: ||log(x1^{-1/2) x2 x1^{-1/2}||_F*/
		virtual realdp DistAF(const Variable &x1, const Variable &x2) const;

		/*Compute the intrinsic representation of a tagnent vector etax, intreta = upper triangle
		of etax. */
		virtual Vector &ObtainIntrEuc(const Variable &x, const Vector &etax, Vector *result) const;

		/*Compute the extrinsic representation of a tagnent vector intretax. Inverse operation of the function ObtainIntrEuc */
		virtual Vector &ObtainExtrEuc(const Variable &x, const Vector &intretax, Vector *result) const;

		/*identity: gf = egf */
		virtual Vector &EucGradToGradEuc(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const;

		/*xix = exix: identity*/
		virtual Vector &EucHvToHvEuc(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const;

		/*compute the distance between two points on the manifold.
		Default: the distance under affine invariant metric: ||log(x1^{-1/2) x2 x1^{-1/2}||_F*/
		virtual realdp DistEuc(const Variable &x1, const Variable &x2) const;
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of SPDMANIFOLD_H */
