
#include "Solvers/ManPG.h"

/*Define the namespace*/
namespace ROPTLIB{

	ManPG::ManPG(const Problem *prob, const Variable *initialx)
	{
		Initialization(prob, initialx);
	};

    ManPG::~ManPG()
    {
    };

	void ManPG::SetDefaultParams()
	{
		SolversNSMPGLS::SetDefaultParams();
        Variant = LSPG_ADALIPSCHITZ;
		SolverName.assign("ManPG");
	};

    void ManPG::Run(void)
    {
        Variable xTemp;
        Vector gfTemp;

        SolversNSMPGLS::Run();

        LSstatus = LSPG_SUCCESS;
        f1 = Prob->f(x1); nf++;
        f2 = f1;
        Prob->Grad(x1, &gf1); ng++;
        ndir0 = 0; ndir1 = 1;
        newslope = 0;
        iter = 0;
        Vector initialstepsize(1 + Max_Iteration), acceptedstepsize(1 + Max_Iteration);
        realdp *timeSeriesptr = timeSeries.ObtainWritePartialData();
        realdp *funSeriesptr = funSeries.ObtainWritePartialData();
        realdp *dirSeriesptr = dirSeries.ObtainWritePartialData();
        realdp *initialstepsizeptr = initialstepsize.ObtainWritePartialData();
        realdp *acceptedstepsizeptr = acceptedstepsize.ObtainWritePartialData();
        
        if (Verbose >= FINALRESULT)
            printf("i:%d,f:%.3e,\n", iter, f1);
            
        if (Verbose >= ITERRESULT)
        {
            timeSeriesptr[iter] = static_cast<realdp > (getTickCount() - starttime) / CLK_PS;
            funSeriesptr[iter] = f1;
            dirSeriesptr[iter] = 0;
            initialstepsizeptr[iter] = 0;
            acceptedstepsizeptr[iter] = 0;
        }
        
        bool isstop = IsStopped();
        /*Start the loop*/
        while ((((! isstop) && iter < Max_Iteration) || iter < Min_Iteration) && LSstatus == LSPG_SUCCESS)
        {
            Mani->TangentSpaceProximalMap(x1, gf1, adavalue, SMtol, SMlambda, Prob, &initD, &SMiter, &SMCGiter, &eta1);
            totalSMiter += SMiter;
            totalSMCGiter += SMCGiter;
            
            ndir1 = sqrt(Mani->Metric(x2, eta1, eta1));
            if(ndir0 == 0)
            { /*initialize ndir0*/
                ndir0 = ndir1;
            }
            
            initialslope = Mani->Metric(x1, gf1, eta1);
            /*initial step size is one*/
            stepsize = 1;
            initiallength = stepsize;
            
            /*Start a line search algorithm.*/
            LinesearchArmijo();
            if(LSstatus == LSPG_MINSTEPSIZE)
            { /*If line search reaches the minimum step size, then solving the proximal mapping more accurate. */
                if(SMlambda == static_cast<realdp> (1e-10) && SMtol == static_cast<realdp> (1e-24))
                {
                    break;
                }
                SMtol = ((SMtol * 0.01 > 1e-24) ? SMtol * 0.01 : 1e-24);
                SMlambda = ((SMlambda * 0.1 > 1e-10) ? SMlambda * 0.1 : 1e-10);
                LSstatus = LSPG_SUCCESS;
                continue;
            }
            
            if(Variant == LSPG_ADALIPSCHITZ)
            {
                if(stepsize == initiallength)
                {
                    adavalue /= 1.01;
                }
                else
                {
                    adavalue = ((adavalue * 1.01 < 1.0) ? adavalue * 1.01 : 1.0);
                }
            }
            
            /*Output debug information if necessary.*/
            if (LSstatus < LSPG_SUCCESS && Verbose > FINALRESULT)
            {
                printf("Linesearch fails! LSstatus:%s\n", LSstatusSetnames[LSstatus].c_str());
            }

            if (std::isnan(f2) || std::isinf(f2)) /*Stop when got a nan or inf*/
            {
                printf("New function value is either nan or inf. Stop!\n");
                break;
            }

            iter++;

            /*Update the Hessian approximation for quasi-Newton methods, etc*/
            UpdateData();
            
            /*Call the function to check whether the stopping criterion is satisfied or not.
            The default function is written in SolversNSM.h and SolversNSM.cpp*/
            isstop = IsStopped();

            if (Verbose >= ITERRESULT)
            {
                /*Output information*/
                if (iter % OutputGap == 0)
                {
                    PrintInfo();
                }
                /*Store debug information in the arrays*/
                timeSeriesptr[iter] = static_cast<realdp>(getTickCount() - starttime) / CLK_PS;
                funSeriesptr[iter] = f2; dirSeriesptr[iter] = ndir1;
                acceptedstepsizeptr[iter] = stepsize;
                initialstepsizeptr[iter] = initiallength;
            }
            /*Switch information at x1 and x2*/
            xTemp = x1; x1 = x2; x2 = xTemp;
            gfTemp = gf1; gf1 = gf2; gf2 = gfTemp;
            
            f1 = f2; /* This need be after the update of pre_funs. It is used in initialstepsize(). */
        }

        ComTime = static_cast<realdp>(getTickCount() - starttime) / CLK_PS;
        if (Verbose >= ITERRESULT)
        {
            lengthSeries = iter + 1;
            SolverInfo.AddToFields("initialstepsize", initialstepsize);
            SolverInfo.AddToFields("acceptedstepsize", acceptedstepsize);
        }
        PrintFinalInfo();
    };

    void ManPG::PrintInfo(void)
    {
        printf("i:%d,f:%.10e,df/f:%.3e,", iter, f2, ((f1 - f2) / std::fabs(f2)));

        printf("|nd|:%.3e,t0:%.2e,t:%.2e,s0:%.2e,s:%.2e,time:%.2e,", ndir1, initiallength, stepsize, initialslope, newslope,  static_cast<realdp>(getTickCount() - starttime) / CLK_PS);

        printf("nf:%d,ng:%d,", nf, ng);
        
        if (nH != 0)
            printf("nH:%d,", nH);
        
        printf("nR:%d,", nR);
        
        if (nV != 0)
            printf("nV(nVp):%d(%d),", nV, nVp);
        
        printf("SMlambda:%.1e,SMtol:%.1e,SMiter:%d,SMCGiter:%d,", SMlambda, SMtol, SMiter, SMCGiter);
        printf("\n");
    };

    void ManPG::PrintFinalInfo(void)
    {
        if (Verbose >= FINALRESULT)
        {
            printf("Iter:%d,f:%.10e,ndir:%.3e,|nd|/|nd0|:%.3e,time:%.2e,nf:%d,ng:%d,nR:%d,", iter, f1, ndir1, ndir1/ndir0, ComTime, nf, ng, nR);
            if (nH != 0)
            {
                printf("nH:%d,", nH);
            }
            if (nV != 0)
            {
                printf("nV(nVp):%d(%d),", nV, nVp);
            }
            printf("SMlambda:%.1e,SMtol:%.1e,totalSMiter:%d,totalSMCGiter:%d,", SMlambda, SMtol, totalSMiter, totalSMCGiter);
            printf("\n");
        }
    };
}; /*end of ROPTLIB namespace*/
