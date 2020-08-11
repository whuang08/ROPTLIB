/*
This file defines the class of the Riemannian Newton method. This code does not follow a particular paper.
The search direction is by applying truncated conjugate gradient to approximately solve the linear
system Hessian[direction] = - gradient.

Solvers --> SolversSM --> SolversSMLS --> RNewton

---- WH
*/

#ifndef RNEWTON_H
#define RNEWTON_H

#include "Solvers/SolversSMLS.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	/* output status of the truncated conjugate gradient. It is an output argument and users don't need to assign this enumerate to any member variable.
	LSSM_NEGCURVTURE: Find negative curvature
	LSSM_LCON: Teminate when the kappa variable takes effect, which indicate the convergence rate is linear.
	LSSM_SCON: Teminate when the theta variable takes effect, which indicate the convergence rate is superlinear.
	LSSM_MAXITER: Teminate when the inner iterations reach the maximum inner iterations specified by the member variable "Max_Inner_Iter"
    LSSM_ERROR: other serious errors
	*/
	enum tCGLSSMstatusSet{ LSSM_NEGCURVTURE, LSSM_LCON, LSSM_SCON, LSSM_MAXITER, LSSM_ERROR, TCGLSSMSTATUSSETLENGTH };

	class RNewton : public SolversSMLS{
	public:
		/*The contructor of RNewton method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.*/
		RNewton(const Problem *prob, const Variable *initialx);

		/*Check whether the parameters about RNewton are legal or not.*/
		virtual void CheckParams(void);

		/*Destructor. Delete temporary Vectors r, z, delta and Hd; Delete the strings of RCGmethods' names*/
		virtual ~RNewton(void);

		/*PARAMSMAP is defined in "def.h" and it is a map from string to realdp, i.e., std::map<std::string, realdp> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);

		/* ===============public parameters below================= */

		/*if useRand is true, then set r to be Hess f(x1)[eta1] + grad f(x1), otherwise, set r to be grad f(x1).
		Default: false*/
		bool useRand;

		/*the maximum iterations allowed for solving the local linear system
		Default: 1000*/
		integer Max_Inner_Iter;

		/*the minimum iterations allowed for solving the local linear system
		Default: 0*/
		integer Min_Inner_Iter;

		/*the theta and kappa are used to check whether stopping criterion of solving local linear system is satisfied or not
		Default: theta 0.1, kappa 0.9*/
		realdp theta;
		realdp kappa;
	protected:
		/*Compute the search direction by using truncated conjugate gradient to approximately solve Hessian[direction] = - gradient*/
		virtual void GetSearchDir(void);

		/*Print information specific to RNewton*/
		virtual void PrintInfo(void);

		/*Run the truncated conjugate gradient method for the local linear system*/
		void tCG_LS(void);
        
        /*Call Solvers::SetProbX function; initialize temporary vectors; and indicate RNewton need action of Hessian.
        INPUT:    prob is the problem which defines the cost function, gradient and possible the action of Hessian
        and specifies the manifold of domain.
        initialx is the initial iterate.*/
        virtual void SetProbX(const Problem *prob, const Variable *initialx);

        /*Setting parameters (member variables) to be default values */
        virtual void SetDefaultParams(void);

		integer innerIter; /*the number of iterations of the truncated conjugate gradient method*/
		tCGLSSMstatusSet tCGLSSMstatus; /*the output status of the truncated conjugate gradient method*/
		std::string *tCGLSSMstatusSetnames; /*the output names of the truncated conjugate gradient method*/
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of RNEWTON_H */
