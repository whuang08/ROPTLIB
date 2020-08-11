/*
This file defines the class for the manifold R_*^{n \times p} / O_p, where R_*^{n \times p} is a n by p full column rank
matrix and O_p is a p-by-p orthogonal group. 

i) This manifold is equivalent to the manifold of the set of symmetric positive semidefinite matricies with rank fixed 
(rank p). The Riemannian metric is
   g_Y(\eta_Y, \xi_Y) = 2 tr(Y^T \eta_Y Y^T \xi_Y + Y^T Y \eta_Y^T \xi_Y)
This Riemannian metric is equivalent to the Euclidean metric on the set of symmetric positive semidefinite matricies with 
rank fixed (rank p)

ii) Riemannian metric: g_x(etax, xix) = \trace(x^T x etax^T xix)
iii) Riemannian metric: \trace(etax^T xix)


Manifold --> SymFixedRankQ

---- WH
*/

#ifndef SYMFIXEDRANKQ_H
#define SYMFIXEDRANKQ_H

#include "Manifolds/Manifold.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	/*Note that not all metrics, retractions and vector transports have been done.*/

	/* Riemannian Metric for the R_*^{n \times p} / O_p manifold:
	EUC: g_x(etax, xix) = 2 \trace(x^T etax x^T xix + x^T x etax^T xix); This metric is equivalent to the Euclidean metric on the symmetric positive semidefinite matrix with rank p
    HGZ, g_x(etax, xix) = \trace(x^T x etax^T xix)
	QEUC, g_x(etax, xix) = \trace(etax^T xix);
	*/
	enum SFRankQMetric { SFRANKQEUC, SFRANKQHGZ, SFRANKQQEUC, SFRANKQMETRICLENGTH };

	/*Vector transport:
	NSOPARALLELIZATION: vector transport by parallelization
	NSOPROJECTION: vector transport by projection
	*/
	enum SFRankQVectorTransport { SFRANKQPARALLELIZATION, SFRANKQPROJECTION, SFRANKQVECTORTRANSPORTLENGTH };

	class SymFixedRankQ : public Manifold{
	public:
		/*Construct the manifold R_*^{n \times p} / O_p, metric: EUC, Vector transport by parallelization*/
		SymFixedRankQ(integer n, integer p = 1);

        /*Default one: Choose Riemannian metric: g_x(etax, xix) = 2 \trace(x^T etax x^T xix + x^T x etax^T xix), vector transport by parallelization */
		void ChooseParamsSet1(void);
        
        /*Choose Riemannian metric: g_x(etax, xix) = \trace(x^T x etax^T xix), vector transport by parallelization */
		void ChooseParamsSet2(void);

		/*Choose Riemannian metric: g_x(etax, xix) = \trace(etax^T xix), vector transport by parallelization */
		void ChooseParamsSet3(void);
        
        /*Choose Riemannian metric: g_x(etax, xix) = 2 \trace(x^T etax x^T xix + x^T x etax^T xix), vector transport by projection */
		void ChooseParamsSet4(void);
        
        /*Choose Riemannian metric: g_x(etax, xix) = \trace(x^T x etax^T xix), vector transport by projection */
		void ChooseParamsSet5(void);

		/*Choose Riemannian metric: g_x(etax, xix) = \trace(etax^T xix), vector transport by projection */
		void ChooseParamsSet6(void);

		virtual ~SymFixedRankQ(void);
        
        /*Randomly generate a point on the manifold.*/
        virtual Variable RandominManifold(void) const;
        
        /*Check whether all the parameters are legal or not.*/
        virtual void CheckParams(void) const;

		/*PARAMSMAP is defined in "def.h" and it is a map from string to realdp, i.e., std::map<std::string, realdp> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);
        
		/*Euclidean metric*/
		virtual realdp Metric(const Variable &x, const Vector &etax, const Vector &xix) const;

		/*Call a member function "IntrProjection" or "ExtrProjection" based on member variable "IsIntrApproach"*/
		virtual Vector &Projection(const Variable &x, const Vector &etax, Vector *result) const;
        
        /*ExtrProjection*/
        virtual Vector &ExtrProjection(const Variable &x, const Vector &etax, Vector *result) const;

        /*ExtrProjection for the QEUC Riemannian metric*/
        virtual Vector &ExtrProjectionQEUC(const Variable &x, const Vector &etax, Vector *result) const;

        /*ExtrProjection for the EUC or HGZ Riemannian metric*/
        virtual Vector &ExtrProjectionEUCorHGZ(const Variable &x, const Vector &etax, Vector *result) const;
        
        /*Householder transformations are used to obtain the intrinsic representation of etax*/
        virtual Vector &ObtainIntr(const Variable &x, const Vector &etax, Vector *result) const;

        /*Householder transformations are used to obtain the extrinsic representation of etax*/
        virtual Vector &ObtainExtr(const Variable &x, const Vector &intretax, Vector *result) const;
        
        /*This is not done yet. */
        virtual Vector &coTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const;
        
        /*The vector transport*/
        virtual Vector &VectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const;
        
        /*The inverse vector transport*/
        virtual Vector &InverseVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result)  const;
        
        /*Compute result = \mathcal{T} * H * \mathcal{T}^{-1}.*/
        virtual LinearOPE &TranHInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, LinearOPE *result) const;
        
		/*For intrinsic representataion of etax. The retraction is the same as the retraction by projection for 
         * the SPSD with fixed rank.*/
		virtual Variable &Retraction(const Variable &x, const Vector &etax, Variable *result) const;

		/*The vector transport by differentiated retraction is the same is the vector tranport by projection for this retraction and manifold.*/
        virtual Vector &DiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/*Euclidean gradient to Riemannian gradient*/
        virtual Vector &EucGradToGrad(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const;

		/*Actioon of Euclidean Hessian to action of Riemannian Hessian*/
        virtual Vector &EucHvToHv(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const;

	protected:

		/*for EUC metric*/
		virtual Vector &EucGradToGradEUC(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const;

		/*for EUC metric*/
		virtual Vector &EucHvToHvEUC(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const;

		/*for EUC metric*/
		virtual Vector &ObtainIntrEUC(const Variable &x, const Vector &etax, Vector *result) const;

		/*for EUC metric*/
		virtual Vector &ObtainExtrEUC(const Variable &x, const Vector &intretax, Vector *result) const;

		/*for QEUC metric*/
		virtual Vector &EucGradToGradQEUC(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const;

		/*for QEUC metric*/
		virtual Vector &EucHvToHvQEUC(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const;

		/*for QEUC metric*/
		virtual Vector &ObtainIntrQEUC(const Variable &x, const Vector &etax, Vector *result) const;

		/*for QEUC metric*/
		virtual Vector &ObtainExtrQEUC(const Variable &x, const Vector &intretax, Vector *result) const;

		/*for HGZ metric*/
		virtual Vector &EucGradToGradHGZ(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const;

		/*for HGZ metric*/
		virtual Vector &EucHvToHvHGZ(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const;

		/*for HGZ metric*/
		virtual Vector &ObtainIntrHGZ(const Variable &x, const Vector &etax, Vector *result) const;

		/*for HGZ metric*/
		virtual Vector &ObtainExtrHGZ(const Variable &x, const Vector &intretax, Vector *result) const;
		
		/*The vector transport by projection*/
		virtual Vector &VectorTransportProj(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const;

		mutable integer n; /*The row*/
		mutable integer p; /*The column*/
		SFRankQMetric metric;
		SFRankQVectorTransport VecTran;
	};
}; /*end of ROPTLIB namespace*/

#endif /* end of SYMFIXEDRANKQ_H */
