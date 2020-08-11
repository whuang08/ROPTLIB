/*
This file defines the class of the Riemannian BFGS method in [HGA2015]
[HGA2018]: Wen Huang, K. A. Gallivan, and P.-A. Absil. A Riemannian BFGS Method without Differentiated Retraction for Nonconvex Optimization Problems
SIAM Journal on Optimization, 28(1):pp.470-495, 2018

Solvers --> SolversSM --> SolversSMLS --> RBFGS

---- WH
*/

#ifndef RBFGS_H
#define RBFGS_H

#include "Solvers/SolversSMLS.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class RBFGS : public SolversSMLS{
	public:
		/*The contructor of RBFGS method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		initialH is the initial inverse Hessian approximation. If the input is nullptr, then the identity is used as the initial approximation.*/
		RBFGS(const Problem *prob, const Variable *initialx, LinearOPE *initialH = nullptr);

		/*Destructor. Delete the vectors and Hessian approximation used in RBFGS, i.e., s and y, H and tildeH*/
		virtual ~RBFGS(void);

		/*Check whether the parameters about RBFGS are legal or not.*/
		virtual void CheckParams(void);

        /*PARAMSMAP is defined in "def.h" and it is a map from string to realdp, i.e., std::map<std::string, realdp> .
        This function is used to set the parameters by the mapping*/
        virtual void SetParams(PARAMSMAP params);

        /*specify whether the cost function is convex or not.
        If yes, then the initial Hessian approximation is a scalar times identity, where the scalar is to
        measure the magnitude of eigenvalues, otherwise, it is identity.
        Default: false*/
        bool isconvex;
        
        /*The same as \epsilon in [LF01, (3.2)]
        [LF01]: D.-H. Li and M. Fukushima. On the global convergence of the BFGS method for nonconvex unconstrained optimization problems.
        SIAM Journal on Optimization, 11(4):1054?064, 2001
        Default: 10^{-4}*/
        realdp nu;
        
        /*The same as \alpha in [LF01, (3.2)]
        [LF01]: D.-H. Li and M. Fukushima. On the global convergence of the BFGS method for nonconvex unconstrained optimization problems.
        SIAM Journal on Optimization, 11(4):1054?064, 2001
        Default: 1*/
        realdp mu;
        
	protected:

		/*Compute the search direction. eta1 = H (-gf1) */
		void GetSearchDir(void);

		/*Update the Hessian approximation if necessary*/
		void UpdateData(void);

		/*Print information specific to RBFGS*/
		virtual void PrintInfo(void);
        
        /*Compute result = H v in RBFGS*/
        virtual Vector &HvRBFGS(const Vector &v, const LinearOPE &H, Vector *result);
        
        /*Update the Hessian approximation for RBFGS if necessary*/
        virtual void UpdateDataRBFGS(void);
        
        /*Initialize the solvers by calling the "SetProbX" and "SetDefultParams" functions.
        INPUT:    prob is the problem which defines the cost function, gradient and possible the action of Hessian
        and specifies the manifold of domain.
        initialx is the initial iterate.
        initialH is the initial inverse Hessian approximation. If the input is nullptr, then the identity is used as the initial approximation.*/
        virtual void Initialization(const Problem *prob, const Variable *initialx, LinearOPE *initialH = nullptr);

        /*Initialize the type of iterates x1, x2 and tangent vectors gf1, gf2, s, y, H and tildeH and obtian the problem and manifold information
        INPUT:    prob is the problem which defines the cost function, gradient and possible the action of Hessian
        and specifies the manifold of domain.
        initialx is the initial iterate.
        initialH is the initial inverse Hessian approximation. If the input is nullptr, then the identity is used as the initial approximation.*/
        virtual void SetProbX(const Problem *prob, const Variable *initialx, LinearOPE *initialH = nullptr);

        /*Setting parameters (member variables) to be default values */
        virtual void SetDefaultParams(void);
        
        bool isupdated; /*Mark whether the (inverse) Hessian approximation is updated*/
        realdp betay, inpsy, inpss, inpyy;  /*betay: \|\xi\| / \|\mathcal{T}_{R_\xi} \xi\| in the locking condition;
                                                phic: the coefficient (1-phic) BFGS + phi DFP in Broyden family method
                                                inpsy: g(s, y); inpss: g(s, s); inpyy: g(y, y); */

        Vector s, y;/*the s, y of current step*/
        LinearOPE H; /*The inverse Hessian approximations for current and next iterations respectively*/
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of RBFGS_H */
