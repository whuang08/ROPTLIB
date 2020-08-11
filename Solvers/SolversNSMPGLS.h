/*
This file defines the common features of proximal gradient based line search methods on manifolds.

Solvers --> SolversNSM --> SolversNSMPGLS

---- WH
*/

#ifndef SOLVERSNSMPGLS_H
#define SOLVERSNSMPGLS_H

#include "Solvers/SolversNSM.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

    enum RPGLSVariant{ LSPG_REGULAR, LSPG_ADALIPSCHITZ, RPGVARIANTLENGTH };

    /* Linesearch status. It is an output argument and users don't need to assign this enumerate to any member variable.
    LSPG_MINSTEPSIZE: line search algorithm reaches the minimum stepsize
    LSPG_SUCCESS: line search algorithm succeeds in finding a point satisfying the line search condition
    */
    enum LSstatusSetPG{ LSPG_MINSTEPSIZE, LSPG_SUCCESS, LSSTATUSSETPGLENGTH };

	class SolversNSMPGLS : public SolversNSM{
	public:
        ~SolversNSMPGLS();
        
		/*Call Solvers::SetProbX function and indicate RCG does not need action of Hessian.
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.*/
		virtual void SetProbX(const Problem *prob, const Variable *initialx);

		/*Setting parameters (member variables) to be default values */
		virtual void SetDefaultParams();
        
        /*Check whether the parameters in RPG are legal or not.*/
        virtual void CheckParams(void);
        
        /*PARAMSMAP is defined in "def.h" and it is a map from string to realdp, i.e., std::map<std::string, realdp> .
        This function is used to set the parameters by the mapping*/
        virtual void SetParams(PARAMSMAP params);
        
        /* Fix or adapt the coefficient in Proximal operator
        Default: LSPG_ADALIPSCHITZ */
        RPGLSVariant Variant;
        
        /*The shrinking parameter in the Armijo-Goldstein condition
        Default: 0.1 */
        realdp LS_ratio;
        
        /*the coefficient of the Armijo condition
        Default: 0.0001 */
        realdp LS_alpha;
        
        /*the minimum stepsize allowed in the linesearch algorithm
        Default: 2e-2*/
        realdp Minstepsize;
	protected:
        /*The Armijo-Goldstein condition by the simple backtracking algorithm*/
        virtual void LinesearchArmijo(void);
        
        /*Evaluate the cost function h(stepsize) = f(R_{x_1}(stepsize * eta1))*/
        virtual realdp h(void);

        /*Evaluate the derivative of cost function h, i.e., h'(stepsize) = \frac{d}{d stepsize} f(R_{x_1}(stepsize * eta1))*/
        virtual realdp dh(void);

        /*When one iteration, some algorithms need to update some information. For example,
            quasi-Newton methods need to update the Hessian approximation. They are done in the following function*/
        virtual void UpdateData(void);
        
        LSstatusSetPG LSstatus;
        Vector eta1;
        Vector eta2;
        
        realdp newslope, initialslope, stepsize, initiallength;
        
        realdp SMlambda;
        realdp SMtol;
        integer SMiter;
        integer totalSMiter;
        integer SMCGiter;
        integer totalSMCGiter;
        integer dimNorVec;
        
        /*for adaptive RPG*/
        realdp adavalue;
        
        /*Data for ProximalMap*/
        Vector initD;
        
        std::string *LSstatusSetnames;    /*This string array is to store the line search status names*/
	};
}; /*end of ROPTLIB namespace*/
#endif /* end of SOLVERSNSMPGLS_H */
