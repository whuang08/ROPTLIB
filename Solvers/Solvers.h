/*
This file defines the abstract base class for all the solvers
It defines the common properties and features of all the solvers

Solvers

---- WH
*/

#ifndef SOLVERS_H
#define SOLVERS_H

#include <iostream>
#include <iomanip>
#include <ctime>
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	/*Specify what information will be output in the algorithm.
	The value should be assigned to the member variable: "Debug",
	NOOUTPUT: no output
	FINALRESULT: final results are outputted
	ITERRESULT: Output information every "OutputGap" iterations, "OutputGap" is a member variable
	DETAILED: Output more than necessary information. Developers can put debug information for this mode.
	The details of output information can be found in Appendix B of the User Manual.
	*/
	enum VERBOSEINFO{ NOOUTPUT, FINALRESULT, ITERRESULT, DETAILED, VERBOSELENGTH };
    extern unsigned long starttime;    /*the start time of running the algorithm*/

	class Solvers{
	public:
		/*Run the algorithm. In this class, this function only initialize debug information and output the name of algorithm.
			This function has been overloaded for all the algorithms*/
		virtual void Run(void);

		/*Check whether the general parameters are legal or not.*/
		virtual void CheckParams(void);

		/*Get the optimizer*/
		inline const Variable GetXopt(void) const { return x1; };

		/*Get the final cost function value*/
		inline realdp Getfinalfun(void) const { return f1; };

		/*Geth the computational wall time of the algorithm*/
		inline realdp GetComTime(void) const { return ComTime; };

		/*Get the number of function evaluations*/
		inline integer Getnf(void) const { return nf; };

		/*Get the number of the gradient evaluations*/
		inline integer Getng(void) const { return ng; };

		/*Get the number of retraction evaluations*/
		inline integer GetnR(void) const { return nR; };

		/*The first time to evaluate the action of vector transport \mathcal{T}_{\eta_x} is usually
		more expensive than the other times to evaluate the action of vector transport \mathcal{T}_{\eta_x}.
		The next two functions are used to
		get the number of action of vector transport (first time)
		get the number of action of vector transport (the other times)
		respectively*/
		inline integer GetnV(void) const { return nV; };
		inline integer GetnVp(void) const { return nVp; };

		/*Get the number of action of Hessian*/
		inline integer GetnH(void) const { return nH; };

		/*Get the number of iterations*/
		inline integer GetIter(void) const { return iter; };

		/*timeSeries, funSeries and gradSeries are three arrays to store the computational time, function value and the norm of gradient
		after each iteration
		The next four functions are used to
		get the number of those three series
		get the computation time series
		get the function values series
		get the dists from soln to the iterates
		respectively. */
		inline integer GetlengthSeries(void) const { return lengthSeries; };
		inline Vector GettimeSeries(void) const { return timeSeries; };
		inline Vector GetfunSeries(void) const { return funSeries; };
		inline Vector GetSolverInfo(void) const { return SolverInfo; };

		/*PARAMSMAP is defined in "def.h" and it is a map from string to realdp, i.e., std::map<std::string, realdp> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);

		/*Beside the three stopping criterion specified by the member variable "Stop_Criterion",
		user also can define a stopping criterion by assigning the following function pointer.
		The code always run this function pointer first if it is not a null pointer. */
		bool(*StopPtr) (const Variable &x, const Vector &funSeries, integer lengthSeries, realdp finalval, realdp initval, const Problem *prob, const Solvers *solver);

		/*Destructor. It is a pure virtual function*/
		virtual ~Solvers(void) = 0;

		/*If a value is less than the Tolerance, then stop the algorithm. The value is specified by the member variable "Stop_Criterion".
		Default: 10^{-6}*/
		realdp Tolerance;

		/*Upper bound for running the algorithm. If excuable time of the algorithm is more the the TimeBound, then the algorithm
		is stopped no matther whether a stopping criterion is satisfied or not.
		Default: 60 * 60 * 24 * 365 (one year);*/
		realdp TimeBound;

		/*Maximum number of iteration. The algorithm will stop if the number of iteration reaches Max_Iteration no matter
		whether a stopping criterion is satisfied or not*/
		integer Max_Iteration;

		/*Minimum number of iteration. The algorithm will keep running if the number of iteration is less than Min_Iteration
		no matther whether a stopping criterion is satisfied or not*/
		integer Min_Iteration;

		/*If the Debug is larger than FINALRESULT, then the information will be outputted every "OutputGap" iterations*/
		integer OutputGap;

		/*Specified the output information of the algorithm. The applicable values are given in the enumerate "DEBUGINFO".*/
		VERBOSEINFO Verbose;

	protected:
		/*Initialize the solvers by calling the "SetProbX" and "SetDefultParams" functions.
			INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
			and specifies the manifold of domain.
			initialx is the initial iterate. */
		virtual void Initialization(const Problem *prob, const Variable *initialx);

		/*Initialize the type of iterates x1, x2 and tangent vectors gf1, gf2 and obtian the problem and manifold information
			INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
			and specifies the manifold of domain.
			initialx is the initial iterate. */
		virtual void SetProbX(const Problem *prob, const Variable *initialx);

		/*Setting parameters (member variables) to be default values */
		virtual void SetDefaultParams(void);

		/*Check whether a stopping criterion is satisfied or not*/
		virtual bool IsStopped(void) = 0;

        /*Print information in every few iterations specific to an algorithm*/
        virtual void PrintInfo(void) = 0;
        
        /*Print last information in an algorithm*/
        virtual void PrintFinalInfo(void) = 0;

		/* algorithm-related variables: */
        Vector gf1, gf2;    /*gf1: gradient at x1, gf2: gradient at x2*/
		Variable x1, x2;	/*x1: current iterate, x2: next iterate*/
		realdp f1, f2;		/*f1: function value at x1, f2: function value at x2*/

		/* Input parameters and functions */
		const Manifold *Mani;	/*The manifold on which the cost function is*/
		const Problem *Prob;	/*The problem which defines the cost function, gradient and probably action of Hessian*/

		/* For debug information */
		integer iter; /*number of iterations*/
		realdp ComTime;	/*the computational time*/
		integer nf, ng, nR, nV, nVp, nH; /*number of function evaluations
										 number of gradient evaluations
										 number of retraction evaluations
										 number of vector transport (See GetnV(void) for details)
										 number of vector transport (See GetnVp(void) for details)
										 number of action of Hessian*/
		Vector timeSeries, funSeries; /*three arrays to store the computational time, function values*/
        integer lengthSeries;        /*the length of above three arrays, i.e., the length of timeSeries, funSeries, distSeries.*/

		Vector SolverInfo;/*Other information for output*/

		std::string SolverName; /*The name of the solver. This is assigned in the constructor function of each derived class*/

		/*new memory for the realdp array Vs, type Vector, with length l*/
		void NewVectors(Vector * &Vs, integer l);

		/*delete memory for the realdp array Vs, type Vector, with length l*/
		void DeleteVectors(Vector * &Vs, integer l);

		/*new memory for the realdp array Xs, type Variable, with length l*/
		void NewVariables(Vector * &Xs, integer l);

		/*delete memory for the realdp array Xs, type Variable, with length l*/
		void DeleteVariables(Vector * &Xs, integer l);
	};

}; /*end of ROPTLIB namespace*/

#endif /* end of SOLVERS_H */
