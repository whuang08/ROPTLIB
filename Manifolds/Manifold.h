/*
This file defines the abstract base class for all the manifolds
It defines the common properties and features of all the manifolds.
Users can write their own manifolds by deriving from this class.
Note that if a function requires some arguments can be the same,
then users need to guarantee that the derived function also
support this property.

Manifold

---- WH
*/

#ifndef MANIFOLD_H
#define MANIFOLD_H

/*#define CHECKMANIFOLDOVERRIDDEN
 */

#include "Problems/Problem.h"
#include "Manifolds/Element.h"
#include <iomanip>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	/*declaration of Problem and ProductManifold. Manifold class should know
	the classes Problem and ProductManifold have been defined somewhere.*/
	class Problem;
	class ProductManifold;

	class Manifold{
	public:
		/*Indicate this class is abstract*/
		virtual ~Manifold(void) = 0;
        
        /*Randomly generate on the unit sphere in the tangent space of the given manifold at a point of the manifold.*/
        virtual Variable RandominManifold(void) const = 0;

		/*Define the Riemannian metric: g_x(etax, xix).
		The default one is the Euclidean metric.*/
		virtual realdp Metric(const Variable &x, const Vector &etax, const Vector &xix) const;
        
        /*Define the action of the linear operator, result = Hx(etax), where Hx is an linear operator on T_x M and
        etax is a tangent vector in T_x M. Note that etax and result must be different, i.e., calling this
        function by
        LinearOPEEta(x, Hx, result, result);
        is illegal.
        The default one is the matrix multiplication, i.e., result = Hx * etax*/
        virtual Vector &LinearOPEEta(const Variable &x, const LinearOPE &Hx, const Vector &etax, Vector *result) const;

		/*etax is in the ambient space. This function projects etax onto the tangent space of x, i.e., result = P_{T_x M} etax
        etax and result can be a same argument, i.e.,
        calling this function by
        Projection(x, result, result);
        is legal.
		Default: result <-- etax*/
		virtual Vector &Projection(const Variable &x, const Vector &etax, Vector *result) const;
        
        /*etax is in the ambient space. This function projects etax onto the tangent space of x, i.e., result = P_{T_x M} etax;
        For this function, both etax and result are represented by intrinsic representations.
        etax and result can be a same argument, i.e.,
        calling this function by
        Projection(x, result, result);
        is legal.
        Default: result <-- etax*/
        virtual Vector &IntrProjection(const Variable &x, const Vector &etax, Vector *result) const;

        /*etax is in the ambient space. This function projects etax onto the tangent space of x, i.e., result = P_{T_x M} etax;
        For this function, both etax and result are represented by extrinsic representations.
        etax and result can be a same argument, i.e.,
        calling this function by
        Projection(x, result, result);
        is legal.
        Default: result <-- etax*/
        virtual Vector &ExtrProjection(const Variable &x, const Vector &etax, Vector *result) const;

        /*Compute result = scalar * etax;etax and result can be a same argument, i.e.,
        calling this function by
        ScaleTimesVector(x, scalar, result, result);
        is legal. */
        virtual Vector &ScalarTimesVector(const Variable &x, const realdp &scalar, const Vector &etax, Vector *result) const;
        
        /*Compute result = scalar * etax + xix; xix and result can be a same argument, i.e.,
        calling this function by
        scalarVectorAddVector(x, scalar, etax, result, result);
        is legal.
        However, scalarVectorAddVector(x, scalar, result, result, result) or scalarVectorAddVector(x, scalar, result, xix, result) is illegal. */
        virtual Vector &ScalarVectorAddVector(const Variable &x, const realdp &scalar, const Vector &etax, const Vector &xix, Vector *result) const;
        
        /*Compute result = scalar1 * etax + scalar2 * xix; xix and result can be a same argument, i.e.,
        calling this function by
        VectorLinearCombination(x, scalar1, etax, scalar2, result, result);
        is legal.
        However, VectorLinearCombination(x, scalar1, result, scalar2, result, result) or VectorLinearCombination(x, scalar1, result, scalar2, xix, result) is illegal. */
        virtual Vector &VectorLinearCombination(const Variable &x, realdp scalar1, const Vector &etax, realdp scalar2, const Vector &xix, Vector *result) const;
        
        /*Compute the retraction result = R_x(etax).
        Default: result = x + etax*/
        virtual Variable &Retraction(const Variable &x, const Vector &etax, Variable *result) const;
        
        /*Compute the inverse of the retraction: result = R_x^{-1}(y).*/
        virtual Vector &InvRetraction(const Variable &x, const Variable &y, Vector *result) const;

		/*Compute the tangent vector result satisfying
			g_y(\mathcal{T}_{R_etax}(xix), xiy) = g_x(xix, result) for all xix \in T_x M,
			where y = R_x(etax), xiy \in T_y M.
			This cotangent vector is used in the RBFGS defined in [RW2012]
			  [RW2012]: W. Ring and B. Wirth. Optimization methods on Riemannian manifolds and their application to shape space.
			   SIAM Journal on Optimization, 22(2):596?27, January 2012
            xiy and result can be a same argument, i.e.,
            calling this function by
            coTangentVector(x, etax, y, result, result);
            is legal.
			Default: result <-- xiy */
		virtual Vector &coTangentVector(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const;

		/*Computes the vector transport by differentiated retraction, i.e., result = \mathcal{T}_{R_etax} (xix)
		if the input IsEtaXiSameDir is true, then this means etax and xix are along a same direction. This implies
		the computations may be simplified.
        xix and result can be a same argument, i.e.,
        calling this function by
        DiffRetraction(x, etax, y, result, result);
        is legal.
		Default: result <-- xix */
		virtual Vector &DiffRetraction(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/*computes beta = \|etax\| / \|\mathcal{T}_{R_etax} etax\|
		Default: beta <-- 1*/
		virtual realdp Beta(const Variable &x, const Vector &etax) const;

		/*compute the distance between two points on the manifold.
		Default: \|x1 - x2\|_F*/
		virtual realdp Dist(const Variable &x1, const Variable &x2) const;

		/*Computes the vector transport, i.e., result = \mathcal{T}_etax (xix)
        xix and result can be a same argument, i.e.,
        calling this function by
        VectorTransport(x, etax, y, result, result);
        is legal.
		Default: result <-- xix*/
		virtual Vector &VectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const;

		/*Computes the inverse vector transport, i.e., result = \mathcal{T}_etax^{-1} (xiy)
        xiy and result can be a same argument, i.e.,
        calling this function by
        InverseVectorTransport(x, etax, y, result, result);
        is legal.
		Default: result <-- xiy*/
		virtual Vector &InverseVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const;

		/*Compute result = H(:, start : end) * \mathcal{T}^{-1}, where H(:, start : end) denotes the matrix formed by columns from "start" to "end".
		The movitation to use "start" and "end" is for the product manifolds. See function TranHInvTran in ProductManifold.h.
        Hx and result can be a same argument, i.e.,
        calling this function by
        HInvTran(x, etax, y, result, start, end, result);
        is legal.
		Default: result <-- Hx*/
		virtual LinearOPE &HInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, integer start, integer end, LinearOPE *result) const;

		/*Compute result = \mathcal{T} * H(start : end, :), where H(start : end, :) denotes the matrix formed by cows from "start" to "end".
		The movitation to use "start" and "end" is for the product manifolds. See function TranHInvTran in ProductManifold.h.
        Hx and result can be a same argument, i.e.,
        calling this function by
        TranH(x, etax, y, result, start, end, result);
        is legal.
		Default: result <-- Hx*/
		virtual LinearOPE &TranH(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, integer start, integer end, LinearOPE *result) const;

		/*Compute result = \mathcal{T} * H * \mathcal{T}^{-1}.
        Hx and result can be a same argument, i.e.,
        calling this function by
        TranHInvTran(x, etax, y, result, result);
        is legal.
		Default: result <-- Hx*/
		virtual LinearOPE &TranHInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, LinearOPE *result) const;

		/*Compute result = Hx + scalar * etax * xix^{\flat}.
        Hx and result can be a same argument, i.e.,
        calling this function by
        HaddScaledRank1OPE(x, result, scalar, etax, xix, result);
        is legal.
		Default: result = Hx + scaler * etax * xix^T*/
		virtual LinearOPE &HaddScaledRank1OPE(const Variable &x, const LinearOPE &Hx, realdp scalar, const Vector &etax, const Vector &xix, LinearOPE *result) const;

		/*Compute etaxflat = etax^{\flat}.
        etax and etaxflat can be a same argument, i.e.,
        calling this function by
        ObtainEtaxFlat(x, etaxflat, etaxflat);
        is legal.
		Default: etaxflat <-- etax*/
		virtual Vector &ObtainEtaxFlat(const Variable &x, const Vector &etax, Vector *result) const;

		/*Compute the intrinsic representation of a tangent vector etax,
		  See details in [Hua2013, Section 9.5]
		  [Hua2013]: W. Huang. Optimization algorithms on Riemannian manifolds with applications.
		  PhD thesis, Florida State University, Department of Mathematics, 2013
        etax and result must be different arguments, i.e.,
        calling this function by
        ObtainIntr(x, result, result);
        is illegal.
		Default: result <-- Zeros*/
		virtual Vector &ObtainIntr(const Variable &x, const Vector &etax, Vector *result) const;

		/*Compute the extrinsic representation of a tangent vector etax,
		   See details in [Hua2013, Section 9.5]
		  [Hua2013]: W. Huang. Optimization algorithms on Riemannian manifolds with applications.
		   PhD thesis, Florida State University, Department of Mathematics, 2013
        etax and result must be different arguments, i.e.,
        calling this function by
        ObtainIntr(x, result, result);
        is illegal.
		Default: result <-- Zeros */
		virtual Vector &ObtainExtr(const Variable &x, const Vector &intretax, Vector *result) const;
        
        /*Compute the intrinsic representation of a vector that is orthogonal to the tangent space at x
        Default: result <-- Zeros */
		virtual Vector &ObtainNorVerIntr(const Variable &x, const Vector &etax, Vector *result) const;
        
        /*Compute the extrinsic representation of a vector that is orthogonal to the tangent space at x
        Default: result <-- Zeros */
		virtual Vector &ObtainNorVerExtr(const Variable &x, const Vector &intretax, Vector *result) const;

        /*Compute result = \argmin_{p \in \T_x M} <p, etax> + 0.5 <p, Weight p> + h(x+p),
        The proximal mapping of h and the jacobi of the proximal mapping are given in "prob" */
        virtual Vector &TangentSpaceProximalMap(Variable &x, const Vector &etax, realdp adavalue, realdp SMtol, realdp SMlambda, const Problem *prob, Vector *inoutinitD, integer *outSMiter, integer *outSMCGiter, Vector *result) const;
        
		/*Get the name of current manifold*/
		inline std::string GetName(void) const { return name; };

		/*Get the dimension of the manifold*/
		inline integer GetIntrDim(void) const { return IntrinsicDim; };

		/*Get the dimension of the ambient space*/
		inline integer GetExtrDim(void) const { return ExtrinsicDim; };

		/*Set the representation of tangent vectors. True if intrinsic is used, false otherwise*/
		inline void SetIsIntrApproach(bool IsIntrAppr) const { IsIntrApproach = IsIntrAppr; };

		/*Get the representation of tangent vectors. True if intrinsic is used, false otherwise*/
		inline bool GetIsIntrinsic(void) const { return IsIntrApproach; };

		/*Get an empty tangent vector using intrinisic representation*/
		inline Vector GetEMPTYINTR(void) const { return EMPTYINTR; };

		/*Get an empty tangent vector using extrinisic representation*/
		inline Vector GetEMPTYEXTR(void) const { return EMPTYEXTR; };
        
        /*Get an empty tangent vector*/
        inline Vector GetEMPTY(void) const { return IsIntrApproach ? EMPTYINTR : EMPTYEXTR; };

		/*Set whether users want to apply the idea in [HGA2015, Section 4.1] to the vector transport defined in the function "VectorTransport".
			[HGA2015]: Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
			SIAM Journal on Optimization, 25(3):1660?685,2015.	.*/
		inline void SetHasHHR(bool inHasHHR) { HasHHR = inHasHHR; };

		/*Check whether the idea in [HGA2015, Section 4.1] is used or not*/
		inline bool GetHasHHR(void) { return HasHHR; };

		/*Check whether all the parameters are legal or not.*/
		virtual void CheckParams(void) const;

		/*Compute the Riemannian gradient from the Euclidean gradient of a function;
		The function is defined in "prob".
		egf is the Euclidean gradient; the output is the Riemannian gradient.
        If the "UseHess" in prob is true, then the egf must be added to x with field name: "EGrad".
        Note that this need be done in all the derived classes of Manifold.
        egf and result can be a same argument, i.e.,
        calling this function by
        EucGradToGrad(x, result, result, prob);
        is legal.
		It is a pure virtual function. It must be overloaded by derived class */
		virtual Vector &EucGradToGrad(const Variable &x, const Vector &egf, const Problem *prob, Vector *result) const = 0;

		/*Compute the Riemannian action of Hessian from the Euclidean action of Hessian of a function;
		The function is defined in "prob".
		the input exix is exix = EucHess [etax], where EucHess is the Euclidean Hessian,
		the output result is result = RieHess[eta], where RieHess is the Riemannian Hessian.
        exix and result can be a same argument, i.e.,
        calling this function by
        EucHvToHv(x, etax, result, result, prob);
        is legal.
		It is a pure virtual function.It must be overloaded by derived class */
		virtual Vector &EucHvToHv(const Variable &x, const Vector &etax, const Vector &exix, const Problem *prob, Vector *result) const = 0;

		/*PARAMSMAP is defined in "def.h" and it is a map from string to realdp, i.e., std::map<std::string, realdp> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);

		/* ===============
		The functions below are used to check the correctness of the Riemannian functions.
		They can be used for people who need to write their own manifold.
		================= */

		/*Check the correctness of the functions: "ObtainIntr" and "ObtainExtr",
		under the assumption that function "ExtrProjection" is correct.
		If the metrics for extrinsic representation and intrinsic representation have been implemented,
		then the outputs of "extr inp" equivalent to "intr inp" means the vector transport is isometric.*/
		virtual void CheckIntrExtr(Variable x) const;

		/*Check the correctness of the Retraction function: "Retraction",
		under the assumption that functions "ExtrProjection", "ObtainIntr", and "ObtainExtr" are correct.*/
		virtual void CheckRetraction(Variable x) const;

		/*Check the correctness of the Retraction function: "DiffRetraction",
		under the assumption that functions "ExtrProjection", "ObtainIntr", "ObtainExtr" and "Retraction" are correct.*/
		virtual void CheckDiffRetraction(Variable x, bool IsEtaXiSameDir = true) const;

		/*Check whether the vector transport and differentiated retraction satisfy the locking condition,
			where the locking condition is define in [HGA2015, (2.8)]
			[HGA2015]: Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
			SIAM Journal on Optimization, 25(3):1660?685, 2015
			under the assumption that functions "ExtrProjection", "ObtainIntr", "ObtainExtr", "Retraction" and "DiffRetraction" are correct.*/
		virtual void CheckLockingCondition(Variable x) const;

		/*Check the correctness of the cotangent vector function: "coTangentVector",
		under the assumption that functions "ExtrProjection", "ObtainIntr", "ObtainExtr", "Retraction" and "DiffRetraction" are correct.*/
		virtual void CheckcoTangentVector(Variable x) const;

		/*Check the isometry of vector transport. This function gives correct result only when extrinsic representation is used.
		If intrinsic representation is used, please use "CheckIntrExtr" instead. */
		virtual void CheckIsometryofVectorTransport(Variable x) const;

		/*Check the isometry of inverse vector transport. This function gives correct result only when extrinsic representation is used.
		If intrinsic representation is used, please use "CheckIntrExtr" instead. */
		virtual void CheckIsometryofInvVectorTransport(Variable x) const;

		/*Check whether \mathcal{T} \circ \mathcal{T}^{-1} is identity or not, where \mathcal{T} is defined by function "VectorTransport"
		and \mathcal{T}^{-1} is defined by function "InverseVectorTransport".*/
		virtual void CheckVecTranComposeInverseVecTran(Variable x) const;

		/*Check whether TranHInvTran is correct or not.*/
		virtual void CheckTranHInvTran(Variable x) const;

		/*Check whether "HaddScaledRank1OPE" is correct or not*/
		virtual void CheckHaddScaledRank1OPE(Variable x) const;

	protected:
		mutable bool HasHHR;/*Mark whether the idea in [HGA2015, Section 4.1] is used or not*/
		mutable bool IsIntrApproach; /*Mark whether intrinsic representation is used for tangent vector or not*/
		std::string name; /*The name of this manifold*/
		integer IntrinsicDim; /*The dimension of this manifold*/
		integer ExtrinsicDim; /*The dimension of the ambient space.*/
		Vector EMPTYINTR; /*The empty tangent vector using intrinsic representation*/
		Vector EMPTYEXTR; /*The empty tangent vector using extrinsic representation*/

		/* The next 5 functions modified above 5 functions such that the locking condition is satisfied.
		The idea in [HGA2015, Section 4.1] is used.
		They are not required to be satisfied.*/

		/*Apply idea in [HGA2015, Section 4.1] to the function "VectorTransport". */
		virtual Vector &LCVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xix, Vector *result) const;

		/*Apply idea in [HGA2015, Section 4.1] to the function "LCInverseVectorTransport". */
		virtual Vector &LCInverseVectorTransport(const Variable &x, const Vector &etax, const Variable &y, const Vector &xiy, Vector *result) const;

		/*Apply idea in [HGA2015, Section 4.1] to the function "LCHInvTran". */
		virtual LinearOPE &LCHInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, integer start, integer end, LinearOPE *result) const;

		/*Apply idea in [HGA2015, Section 4.1] to the function "LCTranH". */
		virtual LinearOPE &LCTranH(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, integer start, integer end, LinearOPE *result) const;

		/*Apply idea in [HGA2015, Section 4.1] to the function "LCTranHInvTran". */
		virtual LinearOPE &LCTranHInvTran(const Variable &x, const Vector &etax, const Variable &y, const LinearOPE &Hx, LinearOPE *result) const;

		/*The function computes unit vectors used in above LC* functions. The idea is in [HGA2015, Section 4.1]. */
		virtual void Obtainnu1nu2forLC(const Variable &x, const Vector &etax, const Variable &y) const;
        
        virtual Vector &ComputeBLambda(const Variable &x, const Vector &d, const Vector &Weight, const Vector &gfx, const Problem *prob, Vector *result) const;

        virtual Vector &EW(const Variable &x, const Vector &BLambda, const Vector &Weight, const Problem *prob, Vector *result) const;

        virtual Vector &GLdW(const Variable &x, const Vector &d, const Vector &Weight, const Vector &BLambda, const Problem *prob, Vector *result) const;

        virtual Vector &calA(const Variable &x, const Vector &DLambda, const Problem *prob, Vector *result) const;

        virtual Vector &calAadj(const Variable &x, const Vector &d, const Problem *prob, Vector *result) const;
                
        Vector &myCG(const Variable &x, const Vector &nb, integer dimNorVec, realdp mu, realdp tau, realdp lambdanFz, integer maxiter, const Vector &Weight, const Vector &BLambda, const Vector &init, const Problem *prob, integer *CGiter, Vector *result) const;
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of MANIFOLD_H */
