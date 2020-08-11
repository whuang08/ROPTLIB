/*
This file defines the class of the limited-memory Riemannian trust-region symmetric rank one update method in [HAG2014]
[HAG2014]: W. Huang, P.-A. Absil, and K. A. Gallivan. A Riemannian symmetric rank-one trustregion method. 
		Mathematical Programming, 150(2):179?16, February 2015.

Solvers --> SolversSM --> SolversSMTR --> LRTRSR1

---- WH
*/

#ifndef LRTRSR1_H
#define LRTRSR1_H

#include "Solvers/SolversSMTR.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class LRTRSR1 : public SolversSMTR{
	public:
		/*The contructor of LRTRSR1 method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.*/
		LRTRSR1(const Problem *prob, const Variable *initialx);

		/*Destructor. Delete the arrays and vectors used in LRTRSR1,
		i.e., vector series S and Y and YMGS, realdp series RHO, matrices PMGQ, SS, and YY, and permutation indices P (see below for details of the series)*/
		virtual ~LRTRSR1(void);

		/*Check whether the parameters about LRBFGS are legal or not.*/
		virtual void CheckParams(void);

		/*Run the algorithm. New memory for vector series S and Y and YMGS, realdp series RHO,
			matrices PMGQ, SS, and YY, and permutation indices P (see below for details of the series). Then call SolversTR::Run*/
		virtual void Run(void);

        /*the number of pairs of s and y in Limited-memory version of quasi-Newton methods;  The same as \ell in [HGA2018]
        [HGA2018]: Wen Huang, K. A. Gallivan, and P.-A. Absil. A Riemannian BFGS Method without Differentiated Retraction for Nonconvex Optimization Problems
        SIAM Journal on Optimization, 28(1):pp.470-495, 2018
        Default: 4*/
        integer LengthSY;
        
        /*If LMrestart is true, then we discard all pairs of (s_k, y_k) in limited-memory quasi-Newton when its number reaches LengthSY;
        Otherwise, we only discard the oldest one and add a new one, therefore the number of pairs of (s_k, y_k) is nondecreasing.*/
        bool LMrestart;
                
	protected:
        /*Call Solvers::SetProbX function and set up the temporary objects for LRTRSR1 algorithm.
        INPUT:    prob is the problem which defines the cost function, gradient and possible the action of Hessian
        and specifies the manifold of domain.
        initialx is the initial iterate.*/
        virtual void SetProbX(const Problem *prob, const Variable *initialx);

        /*Setting parameters (member variables) to be default values */
        virtual void SetDefaultParams(void);
        
        /*PARAMSMAP is defined in "def.h" and it is a map from string to realdp, i.e., std::map<std::string, realdp> .
        This function is used to set the parameters by the mapping*/
        virtual void SetParams(PARAMSMAP params);

        /*Computes the solution to a trust-region subproblem when the quadratic model is defined by
        limited-memory symmetric rank-one (L-SR1) quasi-Newton matrix. The details can be found in [JJR17]
        [JJR17]: On solving L-SR1 trust-region subproblems */
        virtual void tCG_TR(void);

		/*Compute result = H[Eta], where H is the Hessian or the Hessian approximation
			See details in [HAG2014, Section 4]
			[HAG2014]: W. Huang, P.-A. Absil, and K. A. Gallivan. A Riemannian symmetric rank-one trustregion method.
			Mathematical Programming, 150(2):179?16, February 2015.*/
		virtual Vector &HessianEta(const Vector &Eta, Vector *result);

		/*update the pairs of s and y. Add the latest one and remove the oldest one if necessary.*/
		virtual void UpdateData(void);

		/*This function is called when the candidate is accepted. It transports tangent vectors paris of s and y
		from the tangent space at x1 to the tangent space at x2.*/
		virtual void Acceptence(void);

		/*Print information specific to LRTRSR1*/
		virtual void PrintInfo(void);
        
        /*Compute result = H v in LRTRSR1*/
        virtual Vector &HvLRTRSR1(const Vector &Eta, Vector *result);
        
        /*Update the Hessian approximation for LRTRSR1 if necessary*/
        virtual void UpdateDataLRTRSR1(void);

		void ComputeSBySMW(realdp tauStar, const Vector &Psig, const Vector &R, const Vector &PMGQ, const Vector &Psi);
        
        realdp PhiBar_fg(realdp nlambda_min, realdp Delta, const Vector &Lambda, const Vector &a_j, realdp *gf);
        
        realdp PhiBar_fg(realdp nlambda_min, realdp Delta, const Vector &Lambda, const Vector &a_j);

		realdp LocalNewton(realdp x0, integer maxIter, realdp tol, const Vector &Lambda, const Vector &a_j);
        
        bool isupdated; /*Mark whether the (inverse) Hessian approximation is updated*/
        realdp betay, inpsy, inpss, inpyy;  /*betay: \|\xi\| / \|\mathcal{T}_{R_\xi} \xi\| in the locking condition;
                                                phic: the coefficient (1-phic) BFGS + phi DFP in Broyden family method
                                                inpsy: g(s, y); inpss: g(s, s); inpyy: g(y, y); */

        Vector s, y;/*the s, y, and u of current step*/
        
        Vector *S, *Y; /*The stored pairs of s and y*/
        realdp gamma; /*rho: 1/g(s, y) at current iteration, gamma: g(s, y) / g(y, y) for LRBFGS and gamma: g(y, y) / g(s, y) for RTRSR1*/
        integer Currentlength; /*The current length of array S, Y and RHO*/
        integer beginidx; /*The starting index of S, Y and RHO at current iteration*/

        Vector PMGQ; /*PMGQ is the P - gamma Q matrix defined in [HAG2014, (46)]*/
        Vector *YMGS; /*The stored pairs of s and y, and also Y - gamma S*/
        realdp *SS, *SY; /*SS is S^\flat S which is the matrix Q defined in [HAG2014, (46)],
                                SY is the matrix P defined in [HAG2014, (46)] */
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of LRTRSR1_H */
