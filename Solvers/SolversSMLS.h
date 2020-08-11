/*
This file defines the abstract base class for the linesearch-based solvers for smooth objectives
It defines the common properties and features of the linesearch-based solvers for smooth objectives

Solvers --> SolversSM --> SolversSMLS
							
---- WH
*/

#ifndef SOLVERSSMLS_H
#define SOLVERSSMLS_H

#include <iostream>
#include <list>
#include <ctime>
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include "Solvers/SolversSM.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	/* Linesearch algorithms. It should be assigned to the member variable "LineSearch_LS".
	ARMIJO: The Armijo-Goldstein condition. [DS83 Algorithm A6.3.1]
	WOLFE:  The weak Wolfe condition. [DS83 Algorithm A6.3.1mod]
	STRONGWOLFE: The strong Wolfe condition. [NW06 Algorithm 3.5]
	EXACT: The exact line search based on scalar quasi-Newton method
	WOLFELP: The weak Wolfe condition for lipschitz continuous functions
	INPUTFUN: For this option, users can specify their own line search algorithm by assigning the function pointer "LinesearchInput".
	[DS83]: J. E. Dennis and R. B. Schnabel. Numerical methods for unconstrained optimization and nonlinear equations. Springer, New Jersey, 1983
	[NW06]: J. Nocedal and S. J. Wright. Numerical optimization. Springer, second edition, 2006
	*/
	enum LSAlgoSM{ LSSM_ARMIJO, LSSM_WOLFE, LSSM_STRONGWOLFE, LSSM_EXACT, LSSM_INPUTFUN, LSALGOSMLENGTH };

	/* Linesearch status. It is an output argument and users don't need to assign this enumerate to any member variable.
	NOCURVATURE: the second Wolfe condition is not satisfied
	MINSTEPSIZE: line search algorithm reaches the minimum stepsize
	MAXSTEPSIZE: line search algorithm reaches the maximum stepsize
	NONEXACT: exact line search algorithm does not find a point satisfying the inner stopping criterion
	SUCCESS: line search algorithm succeeds in finding a point satisfying the line search condition
	*/
	enum LSstatusSetSM{ LSSM_NOCURVATURE, LSSM_MINSTEPSIZE, LSSM_MAXSTEPSIZE, LSSM_NONEXACT, LSSM_LSERROR, LSSM_SUCCESS, LSSTATUSSETSMLENGTH };

	/*Initial step size in line search algorithm.
	ONESTEP: t0 = one 
	BBSTEP: t0 = g(s, s) / g(s, y), s is the difference of consecutive iterates and y is the difference of the
			gradients at consecutie iterates.
	QUADINT: t0 = [(3.60), NW06]
	QUADINTMOD: t0 = [page 60, NW06]
	[NW06]: J. Nocedal and S. J. Wright. Numerical optimization. Springer, second edition, 2006
	*/
	enum InitStepsizeSetSM{ LSSM_ONESTEP, LSSM_BBSTEP, LSSM_QUADINT, LSSM_QUADINTMOD, LSSM_EXTRBBSTEP, INITSTEPSIZESETSMLENGTH };

	class SolversSMLS : public SolversSM{
	public:
		/*Run the algorithm. This function gives the framework for linesearch based methods*/
		virtual void Run(void);

		/*Check whether the parameters about linesearch algorithms are legal or not.*/
		virtual void CheckParams(void);

		/*PARAMSMAP is defined in "def.h" and it is a map from string to realdp, i.e., std::map<std::string, realdp> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);

		/*Beside the four line search algorithms provided in this library and specified by the member variable "LineSearch_LS",
		user also can define a line search algorithm by assigning the following function pointer.
		User needs to assign LineSearch_LS to be INPUTFUN to call this function. */
		realdp(*LinesearchInput)(integer iter, const Variable &x1, const Vector &exeta1, realdp initialstepsize, realdp initialslope, const Problem *prob, const Solvers *solver);

		/* ===============public parameters below================= */

		/*Line search algorithm. The applicable values are in the enumerate LSAlgo
		Default: ARMIJO */
		LSAlgoSM LineSearch_LS;

		/*If IsPureLSInput is false, then the step size from the user-written line search method, will
		be used as the initial stepsize for the Armijo linesearch. Otherwise, the step size from the user-written line search
		method is used as the final step size.
		Default: false*/
		bool IsPureLSInput;

		/*the coefficient of the Wolfe first condition
		Default: 0.0001 */
		realdp LS_alpha;

		/*the coefficient of the Wolfe second condition
		Default: 0.999*/
		realdp LS_beta;

		/*the minimum stepsize allowed in the linesearch algorithm
		Default: machine eps*/
		realdp Minstepsize;

		/*the maximum stepsize allowed in the linesearch algorithm
		Default: 1000 */
		realdp Maxstepsize;

		/*The coefficient in the Armijo-Goldstein condition
		Default: 0.1 */
		realdp LS_ratio1;

		/*The coefficient in the Armijo-Goldstein condition
		Default: 0.9 */
		realdp LS_ratio2;

		/*Initial stepsize at the first iteration
		Default: 1*/
		realdp Initstepsize;

		/*When the iterate is close to the minimizer (see the annotation of Accuracy), then fixed
		the stepsize to the Finalstepsize if Finalstepsize > 0. If Finalstepsize <= 0, then the proposed
		initial stepsize is used as the accepted stepsize.
		Default: 1*/
		realdp Finalstepsize;
        
        /*the number of previous bb1 stepsize. Used in ABB_min stepsize. See details in [SRTZ2017]
        [SRTZ2017]: On the steplength selection in gradient methods for unconstrained optimization.
        stepsize * id can be used as the initial Hessian approximation in limite-memory quasi-Newton methods
        Default: 0*/
        integer Num_pre_BB;
        
        /*ratio for step size selection. It is the same as \tau in [Algorithm2, SRTZ2017]
        stepsize * id can be used as the initial Hessian approximation in limite-memory quasi-Newton methods
        [SRTZ2017]: On the steplength selection in gradient methods for unconstrained optimization.
        Default: 1, value 0 defines ss/sy stepsize, if the value is 1 and Num_pre_BB is 0, then it defines sy/yy stepsize.*/
        realdp BBratio;

		/* the number of computed functions values. This is used in the nonmonotonic linesearch which
		uses max_{1 \leq i \leq num_pre_funs} (f_{k + 1 - i})
		Default: 1, which defines the standard Armijo line search.*/
		integer Num_pre_funs;

		/* when ngf2 / (ngf0 + Tolerance) <= PreFunsAccuracy, nonmontonic linesearch is used.
		Default: 1e6, which implies the nonmontonic linesearch is not used. */
		realdp PreFunsAccuracy;

		/*Line search algorithm. The applicable values are in the enumerate InitStepsizeSet
		Default: QUADINTMOD (It may alter based on the derived algorithm class) */
		InitStepsizeSetSM InitSteptype;
	protected:
        
		/*Choose what line search algorithm is used*/
		virtual void ChooseLinesearch(void);

		/*Delete objects that are used in this class*/
		virtual ~SolversSMLS(void);

		/*Compute the search direction. It is a pure virtual function.*/
		virtual void GetSearchDir(void) = 0; /* required to be overload in derived class if the derived class is not abstract */

		/*Compute the initial stepsize using [NW06, Page 60]
			[NW06]: J. Nocedal and S. J. Wright. Numerical optimization. Springer, second edition, 2006	*/
		virtual void InitialStepSize(void);
        
        std::list<realdp> pre_BBs; /* Store a few computed BB stepsize ss/sy for initial Hessian approximation using adaptive BB min (ABB_min) idea*/
        
        /*Print information in every few iterations specific to an algorithm*/
        virtual void PrintInfo(void);
        
        /*Print last information in an algorithm*/
        virtual void PrintFinalInfo(void);

		/*Setting parameters (member variables) to be default values */
		virtual void SetDefaultParams(void);

		/*The Armijo-Goldstein condition.[DS83 Algorithm A6.3.1] combined with nonmontone line search
		[DS83] : J.E.Dennis and R.B.Schnabel.Numerical methods for unconstrained optimization and nonlinear equations.Springer, New Jersey, 1983*/
		virtual void LinesearchArmijo(void);

		/*The weak Wolfe condition [DS83 Algorithm A6.3.1mod]
		[DS83] : J.E.Dennis and R.B.Schnabel.Numerical methods for unconstrained optimization and nonlinear equations.Springer, New Jersey, 1983*/
		virtual void LinesearchWolfe(void);

		/*The strong Wolfe condition [NW06 Algorithm 3.5]
		[NW06] : J.Nocedal and S.J.Wright.Numerical optimization.Springer, second edition, 2006*/
		virtual void LinesearchStrongWolfe(void);

		/*Use scalar quasi-Newton method to find a step size as highly accurate as possible*/
		virtual void LinesearchExact(void);

		/*Evaluate the cost function h(stepsize) = f(R_{x_1}(stepsize * eta1))*/
		virtual realdp h(void);

		/*Evaluate the derivative of cost function h, i.e., h'(stepsize) = \frac{d}{d stepsize} f(R_{x_1}(stepsize * eta1))*/
		virtual realdp dh(void);

		/*When one iteration, some algorithms need to update some information. For example,
		quasi-Newton methods need to update the Hessian approximation and nonlinear conjugate gradient
		needs to update the search direction. They are done in the following function*/
		virtual void UpdateData(void);

		/*The "Run" function call this function pointer and the function pointer points to one of the member linesearch functions:
			void LinesearchArmijo(void);
			void LinesearchWolfe(void);
			void LinesearchStrongWolfe(void);
			void LinesearchExact(void); */
		void (SolversSMLS::*Linesearch)(void);

		/* parameters */
		realdp initiallength;	/*The initial stepsize at an iteration*/
		realdp stepsize;		/*The step size*/
		realdp initialslopepre, initialslope, newslope;	/*The slopes for the scalar function h(t) = f(R_{x_1}(t * eta1)) at 0 and at the accepted stepsize*/
		std::list<realdp> pre_funs; /* Store a few computed function values for nonmonotonic line search*/
		LSstatusSetSM LSstatus;	/*The line search status produced by the linesearch algorithm*/
		std::string *LSstatusSetnames;	/*This string array is to store the line search status names*/
	private:
		/*The function used in the strong Wolfe condition. See Algorithm 3.6 in [NW06]
		[NW06] : J.Nocedal and S.J.Wright.Numerical optimization.Springer, second edition, 2006 */
		void Zoom(realdp x1, realdp fx1, realdp slopex1, realdp x2, realdp fx2);
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of SOLVERSSMLS_H */
