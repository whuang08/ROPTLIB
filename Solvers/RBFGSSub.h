/*
This file defines the limited-memory BFGS for locally lipschitz functions on Riemannian manifolds
 
 Solvers --> SolversNSM --> SolversNSMSub --> SolversNSMSubLS --> RBFGSSub

---- WH
*/

#ifndef RBFGSLPSUB_H
#define RBFGSLPSUB_H

#include "Solvers/SolversNSMSubLS.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class RBFGSSub : public SolversNSMSubLS{
	public:

		/*The contructor of RBFGS method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		initialH is the initial inverse Hessian approximation. If the input is nullptr, then the identity is used as the initial approximation.*/
		RBFGSSub(const Problem *prob, const Variable *initialx, LinearOPE *initialH = nullptr);

		/*Destructor. Delete the vectors and Hessian approximation used in RBFGSLPSub, i.e., s and y, H and tildeH*/
		virtual ~RBFGSSub();

		/*Initialize the solvers by calling the "SetProbX" and "SetDefultParams" functions.
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		initialH is the initial inverse Hessian approximation. If the input is nullptr, then the identity is used as the initial approximation.*/
		virtual void Initialization(const Problem *prob, const Variable *initialx, LinearOPE *initialH = nullptr);

		/*Initialize the type of iterates x1, x2 and tangent vectors gf1, gf2, s, y, H and tildeH and obtian the problem and manifold information
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		initialH is the initial inverse Hessian approximation. If the input is nullptr, then the identity is used as the initial approximation.*/
		virtual void SetProbX(const Problem *prob, const Variable *initialx, LinearOPE *initialH = nullptr);

		/*Setting parameters (member variables) to be default values */
		virtual void SetDefaultParams();
        
        /*Compute result = H v in eps-subgradient-based methods*/
        virtual Vector &HvSub(const Vector &v, Vector *result);
        
	protected:
    
        /*Compute result = H v in RBFGS for L-continuous functions*/
        virtual Vector &HvRBFGSSub(const Vector &v, const LinearOPE &H, Vector *result);
    
        /*Update the Hessian approximation for RBFGS for L-continuous functions if necessary*/
        virtual void UpdateDataRBFGSSub(void);

		/*Print information specific to SolversLPSub*/
		virtual void PrintInfo();

		/*Update the Hessian approximation if necessary*/
		void UpdateData(void);
        
        bool isupdated; /*Mark whether the (inverse) Hessian approximation is updated*/
        realdp betay, inpsy, inpss, inpyy;  /*betay: \|\xi\| / \|\mathcal{T}_{R_\xi} \xi\| in the locking condition;
                                                phic: the coefficient (1-phic) BFGS + phi DFP in Broyden family method
                                                inpsy: g(s, y); inpss: g(s, s); inpyy: g(y, y); */

        Vector s, y;/*the s, y of current step*/
        Vector Py; /*the preconditioned y.*/
        LinearOPE H; /*The inverse Hessian approximations for current and next iterations respectively*/
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of RBFGSLPSUB_H */
