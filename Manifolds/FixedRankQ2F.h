/* (2 factor quotient representation for fixed rank manifold)
This file defines the class for the the quotient manifold R_*^{m \times r} \times R_*^{n \times r} / GL(r),
which is equivalent to the fixed-rank manifold R_r^{m \times n},
where R_*^{s \times t} denotes the noncompact Stiefel manifold and GL(r) denotes the generalized linear group.

Riemannian metric is the same as the Euclidean metric in FixedRankE,
g_{(G, H)}( (dG1, dH1), (dG2, dH2) ) = trace( (G dH1^T + dG1 H^T)^T (G dH2^T + dG2 H^T) )
                                     = trace( G^T G dH2^T dH1 + H^T H dG1^T dG2 + H^T dH1 G^T dG2 + dH2^T H dG1^T G )
Note that the above metric is a zero operator in the normal space. One therefore can choose a metric in the normal space
such that any complementary space of the normal space is a horizontal space. Here, we let horizontal space is
H_{G, H} = \{ ([G, G_perp] [K; T], [H, H_perp] [K^T; Q]: K \in R^{r \times r}, T \in R^{(m - r) times r}, Q \in R^{(n - r) \times r}  ) \}


Manifold --> MultiManifolds --> FixedRankQ2F

---- WH
*/

#ifndef FIXEDRANKQ2F_H
#define FIXEDRANKQ2F_H

#include "Manifolds/Euclidean.h"
#include "Manifolds/MultiManifolds.h"

/*Define the namespace*/
namespace ROPTLIB{

	class FixedRankQ2F : public QuotientManifold{
	public:
		/*Construct the quotient manifold R_*^{m \times r} \times R_*^{n \times r} / GL(r),
		which is equivalent to the fixed-rank manifold R_r^{m \times n}.*/
		FixedRankQ2F(integer m, integer n, integer r);
        
        /*Randomly generate a point on the manifold.*/
        virtual Variable RandominManifold(void) const;
        
        /*Check whether all the parameters are legal or not.*/
        virtual void CheckParams(void) const;
        
		/*Delete the manifold.*/
		~FixedRankQ2F(void);
        
        /*Define the Riemannian metric: g_x(etax, xix).
        The default one is the Euclidean metric.*/
        virtual realdp Metric(const Variable &x, const Vector &etax, const Vector &xix) const;
        
        /*etax is in the ambient space. This function projects etax onto the tangent space of x, i.e., result = P_{T_x M} etax */
        virtual Vector &Projection(const Variable &x, const Vector &etax, Vector *result) const;
        
        /*etax is in the total space R_*^{m \times r} \times R_*^{n \times r}. This function projects etax onto the
        tangent space of x. */
        virtual Vector &ExtrProjection(const Variable &x, const Vector &etax, Vector *result) const;

		/*Intrinsic representation has size (m + n - r) r*/
		virtual Vector &ObtainIntr(const Variable &x, const Vector &etax, Vector *result) const;

		/*Compute the extrinsic approach given the intrinsic representation*/
		virtual Vector &ObtainExtr(const Variable &x, const Vector &intretax, Vector *result) const;
        
        /*Compute the tangent vector result satisfying
            g_y(\mathcal{T}_{R_etax}(xix), xiy) = g_x(xix, result) for all xix \in T_x M,
            where y = R_x(etax), xiy \in T_y M. */
        virtual Vector &coTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const;
        
        /*Vector transport by parallelization*/
        virtual Vector &VectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const;
        
        /*Inverse vector transport by parallelization*/
        virtual Vector &InverseVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const;
        
        /*Compute result = \mathcal{T} * H * \mathcal{T}^{-1}.
        Default: result <-- Hx*/
        virtual LinearOPE &TranHInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, LinearOPE *result) const;
        
        /*Compute result = Hx + scalar * etax * xix^{\flat}.
        Default: result = Hx + scaler * etax * xix^T*/
        virtual LinearOPE &HaddScaledRank1OPE(const Variable &x, const LinearOPE &Hx, realdp scalar, const Vector &etax, const Vector &xix, LinearOPE *result) const;

		/*Perform the default retraction of each manifold component*/
		virtual Variable &Retraction(const Variable &x, const Vector &etax, Variable *result) const;

		/*Perform the vector transport by differentiated retraction.*/
		virtual Vector &DiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir = false) const;
		
		/*the Riemannian gradient.*/
		virtual Vector &EucGradToGrad(const Variable &x, const Vector &egf, const Problem *prob, Vector *result)  const;

		/*The action of the Riemannian Hessian.*/
		virtual Vector &EucHvToHv(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const;

	protected:
		mutable integer m; /*the number of row*/
		mutable integer n; /*the number of column*/
		mutable integer r; /*the rank of the matrix*/
        
        void GenerateFieldsX(const Variable &x) const;
        
        /*This is the retraction by projecting (G H^T + dG H^T + G dH^T) onto the set of rank r matrices.
        This is done by projection. Note that this is a second order retraction.*/
        virtual Variable ProjRetraction(const Variable &x, const Vector &etax, Variable *result) const;
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of FIXEDRANKQ2F_H */
