/*
This file defines the abstract base class for all the trust region-based solvers
It defines the common properties and features of all the trust region-based solvers

Solvers --> SolversSM --> SolversSMTR

---- WH
*/

#ifndef SOLVERSSMTR_H
#define SOLVERSSMTR_H

#include "Solvers/SolversSM.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	/* output status of the truncated conjugate gradient. It is an output argument and users don't need to assign this enumerate to any member variable.
	TR_NEGCURVTURE: Find negative curvature
	TR_EXCREGION: The resulting vector is out of the trust region
	TR_LCON: Teminate when the kappa variable takes effect, which indicate the convergence rate is linear.
	TR_SCON: Teminate when the theta variable takes effect, which indicate the convergence rate is superlinear.
	TR_MAXITER: Teminate when the inner iterations reach the maximum inner iterations specified by the member variable "Max_Inner_Iter"
	*/
	enum tCGstatusSetSM{ TRSM_NEGCURVTURE, TRSM_EXCREGION, TRSM_MIN, TRSM_LCON, TRSM_SCON, TRSM_MAXITER, TRSM_ERROR, TCGSTATUSSETSMLENGTH };

	class SolversSMTR : public SolversSM{
	public:
		/*Run the algorithm. This function gives the framework for all the trust region based methods*/
		virtual void Run(void);

		/*Check whether the parameters about trust region algorithms are legal or not.*/
		virtual void CheckParams(void);

		/*PARAMSMAP is defined in "def.h" and it is a map from string to realdp, i.e., std::map<std::string, realdp> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);

		/* ===============public parameters below================= */

		/*if the difference between the local model and the true function is greater than "Acceptence_Rho",
		then accept the candadite iterate. See c in [Step 5 in Algorithm 1, HAG2014].
		[HAG2014]:  W. Huang, P.-A. Absil, and K. A. Gallivan. A Riemannian symmetric rank-one trust-region method. Mathematical Programming.
		Default: 0.1 */
		realdp Acceptence_Rho;

		/*if the local model and the true function does not match well, then the radius is shrinked by "Shrinked_tau",
		See [Step 17 in Algorithm 1, HAG2014]
		Default: 0.25*/
		realdp Shrinked_tau;

		/*if the local model and the true function match good enough, then the radius is magnified by "Magnified_tau",
		See [Step 12 in Algorithm 1, HAG2014]
		Default: 2*/
		realdp Magnified_tau;

		/*Allowed minimum radius of the trust region.
		Default: machine eps*/
		realdp minimum_Delta;

		/*Allowed maximum radius of the trust region.
		Default: 1000*/
		realdp maximum_Delta;

		/*if useRand is true, then set r to be Hess f(x1)[eta1] + grad f(x1), otherwise, set r to be grad f(x1).
		Default: false*/
		bool useRand;

		/*the maximum iterations allowed for solving the local trust region
		Default: 1000*/
		integer Max_Inner_Iter;

		/*the minimum iterations allowed for solving the local trust region
		Default: 0*/
		integer Min_Inner_Iter;

		/*the theta and kappa are used to check whether stopping criterion of solving local model is satisfied or not
		Default values depends on methods. See the User Manual*/
		realdp theta;
		realdp kappa;

		/*The initial radius of the trust region
		Default: 1*/
		realdp initial_Delta;
	protected:
        /*Print information in every few iterations specific to an algorithm*/
        virtual void PrintInfo(void);
        
        /*Print last information in an algorithm*/
        virtual void PrintFinalInfo(void);

		/*Delete objects that are used in this class*/
		virtual ~SolversSMTR(void);

		/*Compute result = H[Eta], where H is the Hessian or the Hessian approximation*/
		virtual Vector &HessianEta(const Vector &Eta, Vector *result) = 0; /* required to be overloaded in derived class */

		/*Compute the initial iterate. The default value is a zero vector.*/
		virtual void InitialVector(void);

		/*Run the truncated conjugate gradient method for the local model*/
		virtual void tCG_TR(void);

		/*Call Solvers::SetProbX function and set up the temporary objects for trust region algorithm.
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		insoln is the true solution. It is not required and only used for research.*/
		virtual void SetProbX(const Problem *prob, const Variable *initialx);

		/*Setting parameters (member variables) to be default values */
		virtual void SetDefaultParams(void);

		/*When one iteration, some algorithms need to update some information. For example,
		quasi-Newton methods need to update the Hessian approximation. They are done in the following function*/
		virtual void UpdateData(void);

		/*This function is called when the candidate is accepted. It can be used to transport tangent vectors and Hessian approximation
		from the tangent space at x1 to the tangent space at x2.*/
		virtual void Acceptence(void);

		/* algorithm-related variables: */
		realdp rho;	/*the difference between the local model and the true function*/
		realdp Delta;	/*the radius of the trust region*/
        Vector Heta2;
		integer innerIter;	/*The number of inner iterations for solving the local model.*/
		tCGstatusSetSM tCGstatusSM; /*The status of solving the local model*/
		std::string *tCGstatusSetSMnames;	/*This string array is to store the trust region status names*/
	};
}; /*end of ROPTLIB namespace*/
#endif
