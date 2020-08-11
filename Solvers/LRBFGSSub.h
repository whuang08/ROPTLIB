/*
This file defines the limited-memory BFGS for locally lipschitz functions on Riemannian manifolds

 Solvers --> SolversNSM --> SolversNSMSub --> SolversNSMSubLS --> LRBFGSSub

---- WH
*/

#ifndef LRBFGSSUB_H
#define LRBFGSSUB_H

#include "Solvers/SolversNSMSubLS.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class LRBFGSSub : public SolversNSMSubLS{
	public:
		/*Run the algorithm. This function gives the framework for the linesearch method*/
		virtual void Run(void);

		/*The contructor of RBFGS method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.*/
		LRBFGSSub(const Problem *prob, const Variable *initialx);

		/*Destructor. Delete the vectors and Hessian approximation used in RBFGSLPSub, i.e., s and y, H and tildeH*/
		virtual ~LRBFGSSub();

		/*Check whether the parameters about RBFGSLPSub are legal or not.*/
		virtual void CheckParams();
        
        /*PARAMSMAP is defined in "def.h" and it is a map from string to realdp, i.e., std::map<std::string, realdp> .
        This function is used to set the parameters by the mapping*/
        virtual void SetParams(PARAMSMAP params);

		/*Initialize the solvers by calling the "SetProbX" and "SetDefultParams" functions.
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.*/
		virtual void Initialization(const Problem *prob, const Variable *initialx);

		/*Initialize the type of iterates x1, x2 and tangent vectors gf1, gf2, s, y, H and tildeH and obtian the problem and manifold information
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.*/
		virtual void SetProbX(const Problem *prob, const Variable *initialx);

		/*Setting parameters (member variables) to be default values */
		virtual void SetDefaultParams();

        /*Compute result = H v in eps-subgradient-based methods*/
        virtual Vector &HvSub(const Vector &v, Vector *result);
        
        /*the number of pairs of s and y in Limited-memory version of quasi-Newton methods;  The same as \ell in [HGA2018]
        [HGA2018]: Wen Huang, K. A. Gallivan, and P.-A. Absil. A Riemannian BFGS Method without Differentiated Retraction for Nonconvex Optimization Problems
        SIAM Journal on Optimization, 28(1):pp.470-495, 2018
        Default: 4*/
        integer LengthSY;
        
	protected:
        
        /*Compute result = H v in LRBFGS for L-continuous functions*/
        virtual Vector &HvLRBFGSSub(const Vector &v, Vector *result);
        
        /*Update the Hessian approximation for LRBFGS for L-continuous functions if necessary*/
        virtual void UpdateDataLRBFGSSub(void);

		/*Print information to LRBFGSSub*/
		virtual void PrintInfo();

		/*Update the Hessian approximation if necessary*/
		void UpdateData(void);
        
        bool isupdated; /*Mark whether the (inverse) Hessian approximation is updated*/
        realdp betay, inpsy, inpss, inpyy;  /*betay: \|\xi\| / \|\mathcal{T}_{R_\xi} \xi\| in the locking condition;
                                                phic: the coefficient (1-phic) BFGS + phi DFP in Broyden family method
                                                inpsy: g(s, y); inpss: g(s, s); inpyy: g(y, y); */
        Vector s, y;/*the s, y, and u of current step*/
        Vector *S, *Y; /*The stored pairs of s and y*/
        realdp *RHO; /*the sequence of 1 / g(s_k, y_k), where Hessian approximation k-th iteration is updated*/
        realdp rho, gamma; /*rho: 1/g(s, y) at current iteration, gamma: g(s, y) / g(y, y) for LRBFGS and gamma: g(y, y) / g(s, y) for RTRSR1*/
        integer Currentlength; /*The current length of array S, Y and RHO*/
        integer beginidx; /*The starting index of S, Y and RHO at current iteration*/
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of LRBFGSSUB_H */
