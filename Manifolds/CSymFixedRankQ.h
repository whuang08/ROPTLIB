/*
This file defines the class for the manifold C_*^{n \times p} / U_p, where C_*^{n \times p} is a n by p full column rank
complex matrix and U_p is a p-by-p unitary group. This manifold is equivalent to the manifold of the set of hermitian positive
semidefinite matricies with rank fixed (rank p).

 Two metrics
 i) This manifold is equivalent to the manifold of the set of Hermitian positive semidefinite matricies with rank fixed
 (rank p). The Riemannian metric is
    g_Y(\eta_Y, \xi_Y) = 2 tr(Y^H \eta_Y Y^H \xi_Y + Y^H Y \eta_Y^H \xi_Y)
 This Riemannian metric is equivalent to the Euclidean metric on the set of Hermitian positive semidefinite matricies with
 rank fixed (rank p)

 ii) Riemannian metric: g_x(etax, xix) = \trace(x^H x etax^H xix)
 Details can be found in [HGZ2015].
	[HGZ2015]:W. Huang, Kyle A. Gallivan, and Xiangxiong Zhang. Solving PhaseLift by low rank Riemannian optimization methods for complex semidefinite constraints.
		U.C.Louvain, UCL-INMA-2015.01, 2015.

Manifold --> CSymFixedRankQ

---- WH
*/

#ifndef CSYMFIXEDRANKQ_H
#define CSYMFIXEDRANKQ_H

#include "Manifolds/Manifold.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

    /*Note that not all metrics, retractions and vector transports have been done.*/

    /* Riemannian Metric for the C_*^{n \times p} / U_p manifold:
    EUC: g_x(etax, xix) = 2 \trace(x^H etax x^H xix + x^H x etax^H xix); This metric is equivalent to the Euclidean metric on the Hermitian positive semidefinite matrix with rank p
    HGZ, g_x(etax, xix) = \trace(x^H x etax^H xix)
    */
    enum CSFRankQMetric { CSFRANKQEUC, CSFRANKQHGZ, CSFRANKQMETRICLENGTH };

    /*Vector transport:
    NSOPARALLELIZATION: vector transport by parallelization
    NSOPROJECTION: vector transport by projection
    */
    enum CSFRankQVectorTransport { CSFRANKQPARALLELIZATION, CSFRANKQPROJECTION, CSFRANKQVECTORTRANSPORTLENGTH };

	class CSymFixedRankQ : public Manifold{
	public:
		/*Construct the manifold C_*^{n \times p} / U_p*/
		CSymFixedRankQ(integer n, integer p = 1);
        
        /*Default one: Choose Riemannian metric: g_x(etax, xix) = 2 \trace(x^H etax x^H xix + x^H x etax^H xix), vector transport by parallelization */
        void ChooseParamsSet1(void);
        
        /*Choose Riemannian metric: g_x(etax, xix) = \trace(x^H x etax^H xix), vector transport by parallelization */
        void ChooseParamsSet2(void);

        /*Choose Riemannian metric: g_x(etax, xix) = 2 \trace(x^H etax x^H xix + x^H x etax^H xix), vector transport by projection */
        void ChooseParamsSet3(void);
        
        /*Choose Riemannian metric: g_x(etax, xix) = \trace(x^H x etax^H xix), vector transport by projection */
        void ChooseParamsSet4(void);

		virtual ~CSymFixedRankQ(void);

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

        /*Householder transformations are used to obtain the intrinsic representation of etax*/
        virtual Vector &ObtainIntr(const Variable &x, const Vector &etax, Vector *result) const;

        /*Householder transformations are used to obtain the extrinsic representation of etax*/
        virtual Vector &ObtainExtr(const Variable &x, const Vector &intretax, Vector *result) const;
        
        /*This is not done yet. */
        virtual Vector &coTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const;
        
        /*The vector transport*/
        virtual Vector &VectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const;
        
        /*The inverse vector transport*/
        virtual Vector &InverseVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const;
        
        /*Compute result = \mathcal{T} * H * \mathcal{T}^{-1}.*/
        virtual LinearOPE &TranHInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, LinearOPE *result) const;
        
        /*For intrinsic representataion of etax. The retraction is the same as the retraction by projection for
         * the SPSD with fixed rank.*/
        virtual Variable &Retraction(const Variable &x, const Vector &etax, Variable *result) const;

        /*The vector transport by differentiated retraction is the same is the vector tranport by projection for this retraction and manifold.*/
        virtual Vector &DiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir = false) const;

        /*Euclidean gradient to Riemannian gradient*/
        virtual Vector &EucGradToGrad(const Variable &x, const Vector &egf, const Problem *prob, Vector *result)  const;

        /*Actioon of Euclidean Hessian to action of Riemannian Hessian*/
        virtual Vector &EucHvToHv(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const;

	protected:
        
        /*for EUC metric*/
        virtual Vector &EucGradToGradEUC(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const;

        /*for EUC metric*/
        virtual Vector &EucHvToHvEUC(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const;
        
        /*for HGZ metric*/
        virtual Vector &EucGradToGradHGZ(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const;

        /*for HGZ metric*/
        virtual Vector &EucHvToHvHGZ(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const;

        /*for EUC metric*/
        virtual Vector &ObtainIntrEUC(const Variable &x, const Vector &etax, Vector *result) const;

        /*for EUC metric*/
        virtual Vector &ObtainExtrEUC(const Variable &x, const Vector &intretax, Vector *result) const;
        
        /*for HGZ metric*/
        virtual Vector &ObtainIntrHGZ(const Variable &x, const Vector &etax, Vector *result) const;

        /*for HGZ metric*/
        virtual Vector &ObtainExtrHGZ(const Variable &x, const Vector &intretax, Vector *result) const;
        
        /*The vector transport by projection*/
        virtual Vector &VectorTransportProj(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const;

        CSFRankQMetric metric;
        CSFRankQVectorTransport VecTran;
		mutable integer n; /*The row*/
		mutable integer p; /*The column*/
	};
}; /*end of ROPTLIB namespace*/

#endif /* end of CSYMFIXEDRANKQ_H */
