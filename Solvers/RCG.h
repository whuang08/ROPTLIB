/*
This file defines the class of the Riemannian nonlinear conjugate gradient method. This code does not follow a particular paper.
It generates all the existing Euclidean nonlinear conjugate gradient methods to the Riemannian setting use generic retraction
and vector transport

Solvers --> SolversSM --> SolversSMLS --> RCG

---- WH
*/

#ifndef RCG_H
#define RCG_H

#include "Solvers/SolversSMLS.h"
#include "Others/def.h"

#undef max

/*Define the namespace*/
namespace ROPTLIB{

	/* Riemannian nonlinear conjugate gradient formulas. It should be assigned to the member variable "RCGmethod".
	The Euclidean formulas can be found in e.g., [NW06, Section 5.2].
	[NW06]: J. Nocedal and S. J. Wright. Numerical optimization. Springer, second edition, 2006
	*/
	enum RCGmethods{ FLETCHER_REEVES, POLAK_RIBIERE_MOD, HESTENES_STIEFEL, FR_PR, DAI_YUAN, HAGER_ZHANG, RCGMETHODSLENGTH };

	class RCG : public SolversSMLS{
	public:
		/*The contructor of RCG method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate. */
		RCG(const Problem *prob, const Variable *initialx);

		/*Destructor. Delete the strings of RCGmethods' names*/
		virtual ~RCG(void);

		/*Check whether the parameters about RCG are legal or not.*/
		virtual void CheckParams(void);

		/*PARAMSMAP is defined in "def.h" and it is a map from string to realdp, i.e., std::map<std::string, realdp> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);

		/* ===============public parameters below================= */

		/*Indicate what formula is used in RCG method*/
		RCGmethods RCGmethod;

	protected:
		/*Compute the search direction based on the RCG forumla.
		Reset the search direction to be negative gradient if the search direction is not sufficiently descent or
        | g( grad f(x_{k+1}), grad f(x_{k}) ) | / g( grad f(x_{k+1}), grad f(x_{k+1}) ) > 0.1 (see [(5.52), NW06]]).*/
		virtual void GetSearchDir(void);

		/*Compute a candadite of the search direction*/
		virtual void UpdateData(void);

		/*Print information specific to RCG*/
		virtual void PrintInfo(void);
        
        /*Call Solvers::SetProbX function and indicate RCG does not need action of Hessian.
        INPUT:    prob is the problem which defines the cost function, gradient and possible the action of Hessian
        and specifies the manifold of domain.
        initialx is the initial iterate.*/
        virtual void SetProbX(const Problem *prob, const Variable *initialx);

        /*Setting parameters (member variables) to be default values */
        virtual void SetDefaultParams(void);

		/*Strings to store the names of formula used in RCG methods*/
		std::string *RCGmethodSetnames;

		/*sigma is the coefficient in - \grad f(x_{k+1}) + sigma eta_k*/
		realdp sigma;

	};
}; /*end of ROPTLIB namespace*/
#endif /* end of RCG_H */
