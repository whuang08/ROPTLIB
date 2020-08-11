/*
This file defines the class of a manifold formed by multiple manifolds.
It can be used to construct a product of manifolds or a quotient manifold whose
total space is a product of manifolds.
It defines the common properties and features of multiple manifolds.
Users can write their own product of manifolds or a quotient manifold by
deriving from this class.

Manifold --> MultiManifolds

---- WH
*/

#ifndef MULTIMANIFOLDS_H
#define MULTIMANIFOLDS_H

/*This is used to help to check whether there is a memory leakage or not.*/
#define CHECKMEMORY

/*ProductManifold and QuotientManifold are just MultiManifolds*/
#define ProductManifold MultiManifolds
#define QuotientManifold MultiManifolds /*The total space of quotient manifold is a product of manifolds*/

#include "Manifolds/Manifold.h"
#include <cstdarg>
#include <map>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class MultiManifolds : public Manifold{
	public:
		/*Constructor of ProductManifold.
		An example of using this constructor to generate St(2,3)^2 \times Euc(2) is:
		integer n = 3, p = 2, m = 2;
		integer numoftypes = 2;
		integer numofmani1 = 2;
		integer numofmani2 = 1;
		Stiefel mani1(n, p);
		Euclidean mani2(m);
		ProductManifold ProdMani(numoftypes, &mani1, numofmani1, &mani2, numofmani2);
		See examples in TestProdStieSumBrockett.cpp and more in test files.	*/
		MultiManifolds(integer numberoftypes, ...);
        
        /*Constructor of ProductManifold
            An example of using this constructor to generate St(2,3)^2 \times Euc(2) is:
            integer n = 3, p = 2, m = 2;
            Stiefel mani1(n, p);
            Euclidean mani2(m);
            Manifold **inmanifolds = new Manifold* [2]; inmanifolds[0] = &mani1; inmanifolds[1] = &mani2;
            integer inpowsinterval = {0, 2, 3};
            ProductManifold ProdMani(inmanifolds, 2, inpowsinterval, 3);
            * The first argument indicates that there are two kinds of manifolds St(2, 3) and Euc(2).
            * The second argement indicates that the length of "inmanifolds" is 2.
            * The third argument indicates the number for each manifold. The number of St(2, 3) is inpowsinterval[1] - inpowsinterval[0] = 2
            and the number of Euc(2) is inpowsinterval[2] - inpowsinterval[1] = 1. Therefore, the product manifold if St(2, 3)^2 \times Euc(2).
            * The fourth argument indicates the number of all manifolds, i.e., 2 Stiefel + 1 Euclidean = 3. */
        MultiManifolds(Manifold **inmanifolds, integer innumoftypes, integer *inpowsinterval);
        
        /*Randomly generate a point on the manifold.*/
        virtual Variable RandominManifold(void) const;

		/*Destructor of ProductManifold. Delete EMPTYINTR, EMPTYEXTR, manifolds, and powsinterval;*/
		virtual ~MultiManifolds(void);

		/*Define the Riemannian metric: g_x(etax, xix).
		The default one is summation of the metrics of all manifolds, i.e.,
		g_x(etax, xix) = Sum g_{x_i}(etax_i, xix_i), where x_i, etax_i, and xix_i are components of x, etax, and xix respectively. */
		virtual realdp Metric(const Variable &x, const Vector &etax, const Vector &xix) const;

		/*Define the action of the linear operator, result = Hx(etax), where Hx is an linear operator on T_x M and
		etax is a tangent vector in T_x M.
		The default one is the matrix multiplication, i.e., result = Hx * etax*/
		virtual Vector &LinearOPEEta(const Variable &x, const LinearOPE &Hx, const Vector &etax, Vector *result) const;
        
        /*Compute result = scalar * etax;*/
        virtual Vector &ScalarTimesVector(const Variable &x, const realdp &scalar, const Vector &etax, Vector *result) const;
        
        /*Compute result = scalar * etax + xix;*/
        virtual Vector &ScalarVectorAddVector(const Variable &x, const realdp &scalar, const Vector &etax, const Vector &xix, Vector *result) const;
        
		/*Compute result = scalar1 * etax + scalar2 * xix; */
        virtual Vector &VectorLinearCombination(const Variable &x, realdp scalar1, const Vector &etax, realdp scalar2, const Vector &xix, Vector *result) const;

		/*etax is in the ambient space. This function projects etax onto the tangent space of x, i.e., result = P_{T_x M} etax
		The component of v and results use the representation specified by each manifold if IsIntrApproach is true.
		Otherwise, the component of v and results use the extrinsic representation.
		The default computation is r_i = P_{T_{x_i} M_i} v_i for all i. */
        virtual Vector &Projection(const Variable &x, const Vector &etax, Vector *result) const;
        
        /*etax is in the ambient space. This function projects etax onto the tangent space of x, i.e., result = P_{T_x M} etax;
        For this function, the components of v and result are represented by extrinsic representations.
        Default function: Let x = (x_1, dots, x_n), v=(v_1, dots, v_n), result = (r_1, dots, r_n).
        The default computation is r_i = P_{T_{x_i} M_i} v_i for all i. */
        virtual Vector &ExtrProjection(const Variable &x, const Vector &etax, Vector *result) const;

		/*Compute the retraction result = R_x(etax).
		Default: Let x = (x_1, dots, x_n), etax=(etax_1, dots, etax_n), result = (r_1, dots, r_n).
		r_i = R_{x_i}(eta_i), for all i. */
        virtual Variable &Retraction(const Variable &x, const Vector &etax, Variable *result) const;
        
        /*Compute the inverse of the retraction: result = R_x^{-1}(y).*/
        virtual Vector &InvRetraction(const Variable &x, const Variable &y, Vector *result) const;

		/*Compute the tangent vector result satisfying
		g_y(\mathcal{T}_{R_etax}(xix), xiy) = g_x(xix, result) for all xix \in T_x M,
		where y = R_x(etax), xiy \in T_y M.
		This cotangent vector is used in the RBFGS defined in [RW2012]
		[RW2012]: W. Ring and B. Wirth. Optimization methods on Riemannian manifolds and their application to shape space.
		SIAM Journal on Optimization, 22(2):596?27, January 2012
		Default: compute the cotangent vector for each component */
        virtual Vector &coTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const;

		/*Computes the vector transport by differentiated retraction, i.e., result = \mathcal{T}_{R_etax} (xix)
		if the input IsEtaXiSameDir is true, then this means etax and xix are along a same direction. This implies
		the computations may be simplied.
		Default: compute the vector transport by differentiated retraction for each component */
        virtual Vector &DiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/*computes beta = \|etax\| / \|\mathcal{T}_{R_etax} etax\|
		Default: etax = (etax_1, dots, eta_n);
		beta = \sqrt(sum \|eta_i\|_2^2) / \sqrt(sum \|\mathcal{T}_{R_{etax_i}} etax_i\|) */
        virtual realdp Beta(const Variable &x, Vector &etax) const;

		/*compute the distance between two points on the manifold.
		Default: sqrt(dist_M1(x1_1, x2_2)^2 + dist_M2(x1_2, x2_2)^2 + ... + dist_M1(x1_n, x2_n)^2)
		where x1_i is the i-th component of the product element.*/
        virtual realdp Dist(const Variable &x1, const Variable &x2) const;

		/*Computes the vector transport, i.e., result = \mathcal{T}_etax (xix)
		Default: compute the vector transport for each component */
        virtual Vector &VectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const;

		/*Computes the inverse vector transport, i.e., result = \mathcal{T}_etax^{-1} (xiy)
		Default: compute the inverse vector transport for each component */
        virtual Vector &InverseVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const;

		/*Compute result = diag(\mathcal{T}_1, dots, \mathcal{T}_n) * H * diag(\mathcal{T}_1^{-1}, dots, \mathcal{T}_n^{-1}).
		Hx and result can be a same argument, i.e.,
		calling this function by
		TranHInvTran(x, etax, y, result, result);
		is legal.
		Default: This is implemented by call the functions "HInvTran" and "TranH" in class "Manifold". */
        virtual LinearOPE &TranHInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, LinearOPE *result) const;

		/*Compute result = Hx + scalar * etax * xix^{\flat}.
		Default: compute xix^{\flat} = (xix_1^\flat, cdots, xix_n^\flat) first, then compute esult = Hx + scalar * etax * xix^{\flat}*/
        virtual LinearOPE &HaddScaledRank1OPE(const Variable &x, const LinearOPE &Hx, realdp scalar, const Vector &etax, const Vector &xix, LinearOPE *result) const;

		/*Each component of etax must use extrinsic representation.
		This function turn the extrinsic representation of a component of etax to intrinsic representation if the
		component manifold has true IsIntrApproach. Otherwise, the component of etax still use extrinsic representation.
		See details in [Hua2013, Section 9.5]
		[Hua2013]: W. Huang. Optimization algorithms on Riemannian manifolds with applications.
		PhD thesis, Florida State University, Department of Mathematics, 2013
		Default: If the "IsIntrApproach" in M_i is true, then call M_i::ObtainIntr(x_i, etax_i, result_i).
		Otherwise, result_i <-- etax_i.*/
        virtual Vector &ObtainIntr(const Variable &x, const Vector &etax, Vector *result) const;

		/*Each component of etax use representation specified by the "IsIntrApproach" of each manifold.
		This function turn the intrinsic representation of a component of etax to extrinsic representation if the
		component manifold has true IsIntrApproach.
		See details in [Hua2013, Section 9.5]
		[Hua2013]: W. Huang. Optimization algorithms on Riemannian manifolds with applications.
		PhD thesis, Florida State University, Department of Mathematics, 2013
		Default: If the "IsIntrApproach" in M_i is true, then call M_i::ObtainExtr(x_i, intretax_i, result_i).
		Otherwise, result_i <-- etax_i.*/
        virtual Vector &ObtainExtr(const Variable &x, const Vector &intretax, Vector *result) const;

		/*Compute the Riemannian gradient from the Euclidean gradient of a function;
		The function is defined in "prob".
		egf is the Euclidean gradient; result is the output which is the Riemannian gradient.
		Default: call the EucGradToGrad for each component manifold. */
        virtual Vector &EucGradToGrad(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const;

		/*Compute the Riemannian action of Hessian from the Euclidean action of Hessian of a function;
		The function is defined in "prob".
		the input exix is exix = EucHess [etax], where EucHess is the Euclidean Hessian,
		the output result is result = RieHess[eta], where RieHess is the Riemannian Hessian.
		Default: call the EucHvToHv for each component manifold. */
        virtual Vector &EucHvToHv(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const;

		/*Check whether all the parameters are legal or not.*/
		virtual void CheckParams(void) const;

		/*Get the idx-th kind of manifolds.*/
		inline Manifold *GetManifold(integer idx) const { if (idx < numoftypes) return manifolds[idx]; return nullptr; };

	protected:
		Manifold **manifolds; /*Store all kinds of manifolds*/
		integer numoftypes; /*The number of kinds of manifolds, i.e., the length of manifolds*/
		integer *powsinterval; /*Each manifold are located in what interval*/
		integer numoftotalmani; /*The number of all manifolds.*/

		/*An example for store St(2, 3) ^ 2 \times Euc(2) is given below:
		(Not exactly C++ code. just give an idea, see an code example in TestProdStieSumBrockett.*)
		Stiefel mani1(3, 2);
		Euclidean mani2(2);
		manifolds = {&mani1, &mani2};  two kinds of manifolds
		numoftypes = 2;  two kinds of manifold
		powsinterval = {0, 2, 3};  from 0 to 2-1 is St(2,3); from 2 to 3-1 is Euc(2);
		numoftotalmani = 3;  (2 Stiefel) + (1 Euclidean) = 3, i.e., 3 manifolds in total.
		*/
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of MULTIMANIFOLDS_H */
