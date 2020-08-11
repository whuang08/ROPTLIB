/*
This file defines the class of a point X on the fixed-rank manifold R_r^{m times n}, which is represented by
an m by n matrix, note that its factors U, D, V such that UDV^T = X, are attached on X with fields names "U",
"D", and "V".
A tangent vector in T_x R_r^{m times n} is also representated by a m by n matrix eta, note that
etax = dot{U} D V^T + U dot{D} V^T + U D dot{V}^T, where dot{U}^T U = 0, dot{V}^T V = 0. Note that dot{U},
dot{D}, and dot{V} are attached on eta with fields names "dU", "dD", and "dV".

Note that the attached fields "U", "D", "V" on a point "X" and "dU", "dD" and "dV" on etax need by modified
accordingly.

The used Riemannian metric is
g(etax, xix) = trace(etax^T xix)
= trace(D^T \dot{U}_1^T \dot{U}_2 D) + \trace(\dot{D}_1^T \dot{D}_2) + \trace(D \dot{V}_1^T \dot{V}_2 D^T),
where etax = dot{U}_1 D V^T + U dot{D}_1 V^T + U D dot{V}_1^T and xix = dot{U}_2 D V^T + U dot{D}_2 V^T + U D dot{V}_2^T.

Manifold --> FixedRankE

---- WH
*/

#ifndef FIXEDRANKE_H
#define FIXEDRANKE_H

#include "Manifolds/Manifold.h"

/*Define the namespace*/
namespace ROPTLIB{

	class FixedRankE : public Manifold{
	public:
		/*Construct the low rank manifold of m by n matrices with rank r.
		It is represented by St(r, m) times R^{r times r} times St(r, n), i.e.,
		X = U D V^T. U \in St(r, m), D \in R^{r times r} and V \in St(r, n).
		Note that D is not necessary a diagonal matrix.*/
		FixedRankE(integer m, integer n, integer r);

		/*Delete the manifold by deleting each component.*/
		~FixedRankE(void);

		/*Riemannian metric*/
        virtual realdp Metric(const Variable &x, const Vector &etax, const Vector &xix) const;
        
        /*Randomly generate a point in the manifold.*/
        virtual Variable RandominManifold() const;

        /*Define the action of the linear operator, result = Hx(etax), where Hx is an linear operator on T_x M and
        etax is a tangent vector in T_x M.
        The default one is the matrix multiplication, i.e., result = Hx * etax*/
        virtual Vector &LinearOPEEta(const Variable &x, const LinearOPE &Hx, const Vector &etax, Vector *result) const;

        /*Compute the intrinsic representation of a vector that is orthogonal to the tangent space at x
        a vector in the normal space is U_perp M V_perp^T, the intrinsic representation is a vector by
        vectorizing the matrix M */
        virtual Vector &ObtainNorVerIntr(const Variable &x, const Vector &etax, Vector *result) const;
        
        /*Compute the extrinsic representation of a vector that is orthogonal to the tangent space at x
         a vector in the normal space is U_perp M V_perp^T, the intrinsic representation is a vector by
         vectorizing the matrix M.*/
         virtual Vector &ObtainNorVerExtr(const Variable &x, const Vector &intretax, Vector *result) const;
        
		/*Perform the default retraction of each manifold component*/
        virtual Variable &Retraction(const Variable &x, const Vector &etax, Variable *result) const;

		/*Compute the tangent vector result satisfying
		g_y(\mathcal{T}_{R_etax}(xix), xiy) = g_x(xix, result) for all xix \in T_x M,
		where y = R_x(etax), xiy \in T_y M.*/
        virtual Vector &coTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const;

		/*Perform the vector transport by differentiated retraction of each individual manifold.*/
        virtual Vector &DiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir = false) const;
		
        /*Compute result = scalar * etax;*/
        virtual Vector &ScalarTimesVector(const Variable &x, const realdp &scalar, const Vector &etax, Vector *result) const;
        
        /*Compute result = scalar * etax + xix;*/
        virtual Vector &ScalarVectorAddVector(const Variable &x, const realdp &scalar, const Vector &etax, const Vector &xix, Vector *result) const;
        
        /*Compute result = scalar1 * etax + scalar2 * xix; If one of etax and xix shares the same memory as result, then it is prefered to let xix be result.*/
        virtual Vector &VectorLinearCombination(const Variable &x, realdp scalar1, const Vector &etax, realdp scalar2, const Vector &xix, Vector *result) const;
        
		/*etax is in the ambient space R^{m \times r} \times R^{r \times r} \times R^{n \times r}. This function projects etax onto the 
		tangent space of x, i.e., result = P_{T_x M} v, where P is based on the selected Riemannian metric*/
        virtual Vector &Projection(const Variable &x, const Vector &etax, Vector *result) const;
        
        /*etax is in the ambient space. This function projects etax onto the tangent space of x, i.e., result = P_{T_x M} etax;
        For this function, both etax and result are represented by extrinsic representations.
        Default: result <-- etax*/
        virtual Vector &ExtrProjection(const Variable &x, const Vector &etax, Vector *result) const;

		/*the Riemannian gradient is obtained by projecting the Euclidean onto the tangent space of x.*/
        virtual Vector &EucGradToGrad(const Variable &x, const Vector &egf, const Problem *prob, Vector *result)  const;

		/*The Riemannian action of the Hessian is obtained by Hess f(x)[etax] = P_x(D grad f(x) [etax]).*/
        virtual Vector &EucHvToHv(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const;

	protected:
		mutable integer m; /*the number of row*/
		mutable integer n; /*the number of column*/
		mutable integer r; /*the rank of the matrix*/
        
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of FIXEDRANKE_H */
