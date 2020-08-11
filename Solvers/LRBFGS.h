/*
This file defines the class of the limited-memory Riemannian BFGS method in [HGA2015]
	[HGA2018]: Wen Huang, K. A. Gallivan, and P.-A. Absil. A Riemannian BFGS Method without Differentiated Retraction for Nonconvex Optimization Problems
    SIAM Journal on Optimization, 28(1):pp.470-495, 2018

 Solvers --> SolversSM --> SolversSMLS --> LRBFGS

---- WH
*/

#ifndef LRBFGS_H
#define LRBFGS_H

#include <cstring>
#include "Solvers/SolversSMLS.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class LRBFGS : public SolversSMLS{
	public:
		/*The contructor of LRBFGS method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		insoln is the true solution. It is not required and only used for research*/
		LRBFGS(const Problem *prob, const Variable *initialx);

		/*Destructor. Delete the arrays and vectors used in LRBFGS, i.e., series S and Y, and series RHO*/
		virtual ~LRBFGS(void);

		/*Check whether the parameters about LRBFGS are legal or not.*/
		virtual void CheckParams(void);

		/*Run the algorithm. New memory for S, Y and RHO. Then call SolversLS::Run*/
		virtual void Run(void);

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
        
        /*the number of pairs of s and y in Limited-memory version of quasi-Newton methods;  The same as \ell in [HGA2018]
        [HGA2018]: Wen Huang, K. A. Gallivan, and P.-A. Absil. A Riemannian BFGS Method without Differentiated Retraction for Nonconvex Optimization Problems
        SIAM Journal on Optimization, 28(1):pp.470-495, 2018
        Default: 4*/
        integer LengthSY;
        
        /*If LMrestart is true, then we discard all pairs of (s_k, y_k) in limited-memory quasi-Newton when its number reaches LengthSY;
        Otherwise, we only discard the oldest one and add a new one, therefore the number of pairs of (s_k, y_k) is nondecreasing.*/
        bool LMrestart;
                
	protected:

		/*Compute the search direction. [HGA2018]
			[HGA2018]: Wen Huang, K. A. Gallivan, and P.-A. Absil. A Riemannian BFGS Method without Differentiated Retraction for Nonconvex Optimization Problems
            SIAM Journal on Optimization, 28(1):pp.470-495, 2018
			*/
		virtual void GetSearchDir(void);

		/*update the pairs of s and y. Add the latest one and remove the oldest one if necessary.
		transport them to the tangent space of x2*/
		virtual void UpdateData(void);

		/*Print information specific to LRBFGS*/
		virtual void PrintInfo(void);
        
        /*Compute result = H v in LRBFGS*/
        virtual Vector &HvLRBFGS(const Vector &v, Vector *result);
        
        /*initial Hessian approximation in limited-memory BFGS methods. It is a scalar times identity.*/
        virtual realdp InitialHessian(realdp inpss, realdp inpsy, realdp inpyy);
        
        /*Update the Hessian approximation for LRBFGS if necessary*/
        virtual void UpdateDataLRBFGS(void);
        
        /*Call Solvers::SetProbX function and set up the temporary objects for LRBFGS algorithm.
        INPUT:    prob is the problem which defines the cost function, gradient and possible the action of Hessian
        and specifies the manifold of domain.
        initialx is the initial iterate.*/
        virtual void SetProbX(const Problem *prob, const Variable *initialx);

        /*Setting parameters (member variables) to be default values */
        virtual void SetDefaultParams(void);
        
        bool isupdated; /*Mark whether the (inverse) Hessian approximation is updated*/
        realdp betay, inpsy, inpss, inpyy;  /*betay: \|\xi\| / \|\mathcal{T}_{R_\xi} \xi\| in the locking condition;
                                                phic: the coefficient (1-phic) BFGS + phi DFP in Broyden family method
                                                inpsy: g(s, y); inpss: g(s, s); inpyy: g(y, y); */
        Vector s, y;/*the s, y, and u of current step*/
        Vector Py; /*the preconditioned y.*/
        Vector *S, *Y; /*The stored pairs of s and y*/
        realdp *RHO; /*the sequence of 1 / g(s_k, y_k), where Hessian approximation k-th iteration is updated*/
        realdp rho, gamma; /*rho: 1/g(s, y) at current iteration, gamma: g(s, y) / g(y, y) for LRBFGS and gamma: g(y, y) / g(s, y) for RTRSR1*/
        integer Currentlength; /*The current length of array S, Y and RHO*/
        integer beginidx; /*The starting index of S, Y and RHO at current iteration*/
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of RBROYDENFAMILY_H */
