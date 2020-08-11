/*
This file defines the abstract base class for all the solvers for smooth objectives
It defines the common properties and features of all the solvers for smooth objectives

Solvers --> SolversSM

---- WH
*/

#ifndef SOLVERSSM_H
#define SOLVERSSM_H

#include <iostream>
#include <iomanip>
#include <ctime>
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include "Solvers/Solvers.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	/*The algorithm is stopped when a value (specified by ther parameter) is less than the "Tolerance" (a member variable)
	The value should be assigned to the member variable: "Stop_Criterion" and the applicable values are
	SM_FUN_REL: |f_k - f_{k+1}| / max(|f_k|, 1)
	SM_GRAD_F: \|gf_k\|
	SM_GRAD_F_0: \|gf_k\| / \|gf_0\|*/
	enum StopCritSM{ SM_FUN_REL, SM_GRAD_F, SM_GRAD_F_0, STOPCRITSMLENGTH };

	class SolversSM : public Solvers{
	public:
		/*Run the algorithm.*/
		virtual void Run(void);

		/*Check whether the general parameters are legal or not.*/
		virtual void CheckParams(void);

		/*Get the norm of the final gradient*/
		inline realdp Getnormgf(void) const { return ngf1; };

		/*Get the norm of the gradient at final iterate over the norm of the gradient at initiate iterate*/
		inline realdp Getnormgfgf0(void) const { return ngf1 / ngf0; };

		/*gradSeries is an array to store the norm of gradient after each iteration
		get the norm of gradient series respectively. */
		inline Vector GetgradSeries(void) const { return gradSeries; };

		/*PARAMSMAP is defined in "def.h" and it is a map from string to realdp, i.e., std::map<std::string, realdp> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);

		/*Destructor. It is a pure virtual function*/
		virtual ~SolversSM(void) = 0;

		/*member variable: specify the stopping criterion. The applicable values are defined in enumerate "StopCrit"
		Default: GRAD_F_0*/
		StopCritSM Stop_Criterion;

		/*If the norm of current gradient over the norm of the initial gradient is less than Accuracy,
		then i) in line search strategy: the step size is fixed to be the member variable "Finalstepsize";
		 ii) in trust-region strategy: the candidate is always accepted if the minimizer in the local model is
		 in the constrained ball.
		Defaut: 0*/
		realdp Accuracy;

	protected:
		/*Initialize the solvers by calling the "SetProbX" and "SetDefultParams" functions.
			INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
			and specifies the manifold of domain.
			initialx is the initial iterate. */
		virtual void Initialization(const Problem *prob, const Variable *initialx);

		/*Setting parameters (member variables) to be default values */
		virtual void SetDefaultParams(void);
        
        /*Initialize the type of iterates x1, x2 and tangent vectors gf1, gf2 and obtian the problem and manifold information
            INPUT:    prob is the problem which defines the cost function, gradient and possible the action of Hessian
            and specifies the manifold of domain.
            initialx is the initial iterate. */
        virtual void SetProbX(const Problem *prob, const Variable *initialx);

		/*When one iteration, some algorithms need to update some information. For example,
			quasi-Newton methods need to update the Hessian approximation and nonlinear conjugate gradient
			needs to update the search direction. They are done in the following function*/
		virtual void UpdateData(void) = 0;

		/*Check whether a stopping criterion is satisfied or not*/
		virtual bool IsStopped(void);
        
        /*Print information in every few iterations specific to an algorithm*/
        virtual void PrintInfo(void);
        
        /*Print last information in an algorithm*/
        virtual void PrintFinalInfo(void);

		/* algorithm-related variables: */
		Vector Pgf1, Pgf2;	/*Pgf1: preconditioned gradient at x1, Pgf2: preconditioned gradient at x2, used in RCG and LRBFGS, LRTRSR1*/
		realdp ngf0, ngf1, ngf2; /*ngf0: the norm of gradient at the initial iterate, ngf1: the norm of the gradient at x1, ngf2: the norm of the gradient at x2*/
		Vector eta1, eta2; /*In Line search-based methods, eta1 is the search direction. eta2 is stepsize * eta1. In trust region-based methods, eta1 is the
                             initial guess for solving the local model. eta2 is the result produced by solving the local model.*/

		/* For debug information */
		Vector gradSeries; /*an array to store the norm of gradient after each iteration */
	};

}; /*end of ROPTLIB namespace*/

#endif /* end of SOLVERSSM_H */
