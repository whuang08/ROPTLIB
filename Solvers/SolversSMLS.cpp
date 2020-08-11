
#include "Solvers/SolversSMLS.h"

/*Define the namespace*/
namespace ROPTLIB{

	void SolversSMLS::Run(void)
	{
        bool nonmonoticLS = false;
        Variable xTemp(x1);
        Vector gfTemp = Prob->GetDomain()->GetEMPTY();
		pre_funs.clear();

		SolversSM::Run();

		ChooseLinesearch();

		LSstatus = LSSM_SUCCESS;
		f1 = Prob->f(x1); nf++;
		f2 = f1;
		Prob->Grad(x1, &gf1); ng++;
		ngf0 = sqrt(Mani->Metric(x1, gf1, gf1));
		ngf1 = ngf0; ngf2 = ngf1;
		newslope = 0;
        initialslopepre = 0;
		iter = 0;
        Vector initialstepsize(1 + Max_Iteration), acceptedstepsize(1 + Max_Iteration);
        realdp *timeSeriesptr = timeSeries.ObtainWritePartialData();
        realdp *funSeriesptr = funSeries.ObtainWritePartialData();
        realdp *gradSeriesptr = gradSeries.ObtainWritePartialData();
        realdp *initialstepsizeptr = initialstepsize.ObtainWritePartialData();
        realdp *acceptedstepsizeptr = acceptedstepsize.ObtainWritePartialData();
        
        if (Verbose >= FINALRESULT)
            printf("i:%d,f:%.3e,|gf|:%.3e,\n", iter, f1, ngf1);
		if (Verbose >= ITERRESULT)
		{
			timeSeriesptr[iter] = static_cast<realdp>(getTickCount() - starttime) / CLK_PS;
			funSeriesptr[iter] = f1;
			gradSeriesptr[iter] = ngf1;
			initialstepsizeptr[iter] = 0;
			acceptedstepsizeptr[iter] = 0;
		}
		bool isstop = IsStopped();
        Vector exeta1(Mani->GetEMPTYEXTR());
		/*Start the loop*/
		while ((((! isstop) && iter < Max_Iteration) || iter < Min_Iteration) && LSstatus == LSSM_SUCCESS)
		{
			GetSearchDir(); /* Obtain search direction eta1 */
			initialslope = Mani->Metric(x1, gf1, eta1);
			/*Compute initial step size for the next iteration*/
			InitialStepSize();
			initiallength = stepsize;
			if (Verbose >= ITERRESULT)
				initialstepsizeptr[iter + 1] = initiallength;
            
			/*If accurate enough, then a fixed stepsize is chosen.*/
			if (ngf1 / (ngf0 + Tolerance) < Accuracy)
			{
                stepsize = (Finalstepsize > 0) ? Finalstepsize : initiallength;
				f2 = h();
				Prob->Grad(x2, &gf2); ng++;
			}
			else
			{
				/* Call the specified linesearch algorithm.
				Note that in the linesearch algorithm, we need to obtain
				accepted stepsize, eta2=stepsize*eta1, x2 = R_{x_1}(eta_2), f2 = f(x2), and gf2 = grad f(x_2) */
				if (LineSearch_LS == LSSM_INPUTFUN)
				{
                    stepsize = (Prob->GetDomain()->GetIsIntrinsic()) ?
                        LinesearchInput(iter, x1, Prob->GetDomain()->ObtainExtr(x1, eta1, &exeta1), stepsize, initialslope, Prob, this) :
                        LinesearchInput(iter, x1, eta1, stepsize, initialslope, Prob, this);
                    
					stepsize = (stepsize < std::numeric_limits<realdp>::epsilon()) ? initiallength : stepsize;
					initiallength = stepsize;
					if (!IsPureLSInput)
					{
						SolversSMLS::LinesearchArmijo();
					}
					else
					{
						f2 = h();
						Prob->Grad(x2, &gf2); ng++;
					}
				}
				else
				{
					(this->*Linesearch)();
				}
			}
            
			/*Output debug information if necessary.*/
			if (LSstatus < LSSM_SUCCESS && Verbose > FINALRESULT)
			{
				printf("Linesearch fails! LSstatus:%s\n", LSstatusSetnames[LSstatus].c_str());
			}

			if (std::isnan(f2) || std::isinf(f2)) /*Stop when got a nan or inf*/
			{
				printf("New function value is either nan or inf. Stop!\n");
				break;
			}
			iter++;

			/*Update the Hessian approximation for quasi-Newton methods 
				or obtain search direction candadite for Riemannian nonlinear conjugate gradient*/
			UpdateData();

			/*norm of the gradient at x2*/
			ngf2 = sqrt(Mani->Metric(x2, gf2, gf2));
            
            /*Call the function to check whether the stopping criterion is satisfied or not.
            The default function is written in Solvers.h and Solvers.cpp*/
            isstop = IsStopped();

			if (Verbose >= ITERRESULT)
			{
				/*Output information*/
				if (iter % OutputGap == 0)
					PrintInfo();
                
				/*Store debug information in the arrays*/
				timeSeriesptr[iter] = static_cast<realdp>(getTickCount() - starttime) / CLK_PS;
				funSeriesptr[iter] = f2; gradSeriesptr[iter] = ngf2;
				acceptedstepsizeptr[iter] = stepsize;
			}
			/*Switch information at x1 and x2*/
            xTemp = x1;
            x1 = x2;
            x2 = xTemp;
            gfTemp = gf1; gf1 = gf2; gf2 = gfTemp;
            
            if(ngf2 / (ngf0 + Tolerance) <= PreFunsAccuracy)
                nonmonoticLS = true;
            
			pre_funs.push_front(f1);
			if ((pre_funs.size() > Num_pre_funs || !nonmonoticLS) && pre_funs.size() > 1)
				pre_funs.pop_back();
			f1 = f2; /* This need be after the update of pre_funs. It is used in initialstepsize(). */
			ngf1 = ngf2;
            initialslopepre = initialslope;
		}

		ComTime = static_cast<realdp>(getTickCount() - starttime) / CLK_PS;
		if (Verbose >= ITERRESULT)
		{
			lengthSeries = iter + 1;
            SolverInfo.AddToFields("initialstepsize", initialstepsize);
            SolverInfo.AddToFields("acceptedstepsize", acceptedstepsize);
		}
        if (Verbose >= FINALRESULT)
            PrintFinalInfo();
	};
	
	void SolversSMLS::ChooseLinesearch(void)
	{
		/*If the linesearch condition Armijo-Glodstein condition, then the locking conidition is not necessary and don't output warning.
		If the pair of retraction and vector transport satisfies the locking condition, then output warning.
		If the idea in [Section 4.1, HGA2015] is used, then the locking condition is satisfied and don't output the warning.
		[HGA2015]: Wen Huang, K. A. Gallivan, and P.-A. Absil, A Broyden Class of Quais-Newton Methods for Riemannian Optimization,
		SIAM on Journal Optimization, (25)3, 1660-1685, 2015 */
		if (LineSearch_LS != LSSM_ARMIJO && !Prob->GetDomain()->GetHasHHR() && Verbose >= ITERRESULT)
			printf("Line search algorithm requires the locking condition.\n");

		/*Choose the linesearch algorithm used in the algorithm*/
		if (LineSearch_LS == LSSM_ARMIJO)
			Linesearch = &SolversSMLS::LinesearchArmijo;
		else
		if (LineSearch_LS == LSSM_WOLFE)
			Linesearch = &SolversSMLS::LinesearchWolfe;
		else
		if (LineSearch_LS == LSSM_STRONGWOLFE)
			Linesearch = &SolversSMLS::LinesearchStrongWolfe;
		else
		if (LineSearch_LS == LSSM_EXACT)
			Linesearch = &SolversSMLS::LinesearchExact;
		else
		if (LineSearch_LS == LSSM_INPUTFUN)
		{
			if (LinesearchInput == nullptr)
			{
				printf("Error: linesearch function pointer does not exist!\n");
				return;
			}
		}
		else
		{
			/*If the line search algorithm is not specified, then use the Armijo-Goldstein condition.*/
			if (Verbose >= FINALRESULT)
			{
				printf("Warning: linesearch algorithm does not exist!\n");
				printf("Use linesearch algorithm with Armijo-Goldstein conditions!\n");
			}
			Linesearch = &SolversSMLS::LinesearchArmijo;
		}
	};

	void SolversSMLS::InitialStepSize(void)
	{
        Vector s = Prob->GetDomain()->GetEMPTY(), y(s);
		realdp BB1 = 0, BB2 = 0, quadintstepsize=0;
		if (iter == 0)
			stepsize = Initstepsize;
		else
		{
            quadintstepsize = static_cast<realdp> (2) * (f1 - pre_funs.front()) / initialslopepre;
			switch (InitSteptype)
			{
			case LSSM_QUADINTMOD:
				stepsize = static_cast<realdp> (1.01 * 2.0) * (f1 - pre_funs.front()) / initialslopepre;
				stepsize = (stepsize > 1) ? 1 : stepsize;
				break;
			case LSSM_BBSTEP:
				/*Since x1 and x2 are swapped, gf1 and gf2 are swapped, we have the following formula.*/
                    
                Mani->VectorTransport(x2, eta2, x1, eta2, &s);
                Mani->VectorTransport(x2, eta2, x1, gf2, &y);
                Mani->VectorLinearCombination(x1, 1, gf1, -1, y, &y);
				BB1 = Mani->Metric(x2, s, s) / Mani->Metric(x2, s, y);
				BB2 = Mani->Metric(x2, s, y) / Mani->Metric(x2, y, y);
				pre_BBs.push_front((BB2 < std::numeric_limits<realdp>::epsilon()) ? Initstepsize / ngf1 : BB2);
				if (pre_BBs.size() > Num_pre_BB)
					pre_BBs.pop_back();

				if(BBratio == 0)
					stepsize = BB1;
				else
				if (BBratio == 1 && Num_pre_BB == 0)
				{
					stepsize = BB2;
				}
				else
				{
					realdp minpreBB = BB2;
					std::list<realdp>::iterator j = pre_BBs.begin();
					for (integer i = 0; i < Num_pre_BB && j != pre_BBs.end(); i++, j++)
					{
						if (minpreBB > *j)
							minpreBB = *j;
					}
					if (minpreBB / BB1 < BBratio)
						stepsize = minpreBB;
					else
						stepsize = BB1;
				}
				break;
			case LSSM_EXTRBBSTEP:
				/*Since x1 and x2 are swapped, we have the following formula.*/
                Mani->VectorLinearCombination(x1, 1, x1, -1, x2, &s);
				if (Mani->GetIsIntrinsic())
				{
                    Vector gf1extr(Mani->GetEMPTYEXTR()), gf2extr(Mani->GetEMPTYEXTR());
					Mani->ObtainExtr(x1, gf1, &gf1extr);
					Mani->ObtainExtr(x2, gf2, &gf2extr);
                    Mani->VectorLinearCombination(x1, 1, gf1extr, -1, gf2extr, &y);
				}
				else
				{
                    Mani->VectorLinearCombination(x1, 1, gf1, -1, gf2, &y);
				}
				BB1 = Mani->Metric(x2, s, s) / Mani->Metric(x2, s, y);
				BB2 = Mani->Metric(x2, s, y) / Mani->Metric(x2, y, y);

				pre_BBs.push_front((BB2 < std::numeric_limits<realdp>::epsilon()) ? Initstepsize / ngf1 : BB2);
				if (pre_BBs.size() > Num_pre_BB)
					pre_BBs.pop_back();

				if(BBratio == 0)
					stepsize = BB1;
				else
				if (BBratio == 1 && Num_pre_BB == 0)
				{
					stepsize = BB2;
				}
				else
				{
					realdp minpreBB = BB2;
					std::list<realdp>::iterator j = pre_BBs.begin();
					for (integer i = 0; i < Num_pre_BB && j != pre_BBs.end(); i++, j++)
					{
						if (minpreBB > *j)
							minpreBB = *j;
					}
					if (minpreBB / BB1 < BBratio)
						stepsize = minpreBB;
					else
						stepsize = BB1;
				}
				break;
			case LSSM_ONESTEP:
				stepsize = 1;
				break;
			case LSSM_QUADINT:
                stepsize = quadintstepsize;
				break;
			default:
				printf("InitSteptype is incorrect. Use one instead.\n");
				stepsize = 1;
			};
			/*Safeguard for the initial stepsize*/
            stepsize = (stepsize < quadintstepsize / 100) ? quadintstepsize : stepsize;
			stepsize = (stepsize < std::numeric_limits<realdp>::epsilon()) ? Initstepsize / ngf1 : stepsize;
		}
	};

	void SolversSMLS::LinesearchArmijo(void)
	{
		LSstatus = LSSM_SUCCESS;
		f2 = h();
		realdp maxpref = f1; 
		std::list<realdp>::iterator j;
		std::list<realdp>::iterator jend;
		if (pre_funs.size() > 0)
		{
			j = pre_funs.begin();
			jend = pre_funs.end();
			jend--;
		}
		for (integer i = 0; i < Num_pre_funs && j != jend; i++, j++)
		{
			if (maxpref < *j)
				maxpref = *j;
		}
		/* simple backtracking */
		if (LS_ratio2 <= LS_ratio1)
		{
			realdp LS_ratio = LS_ratio1;
			while (maxpref - f2 < -LS_alpha * initialslope * stepsize)
			{
				stepsize *= LS_ratio;
				if (stepsize < Minstepsize)
				{
					if (Verbose > FINALRESULT)
					{
						printf("Warning: step size reaches the minimum: %3.2e!\n", Minstepsize);
					}
					LSstatus = LSSM_MINSTEPSIZE;
					break;
				}
				f2 = h();
			}
            newslope = dh();
			return;
		}

		/* use polynomial interplation */
		realdp prestepsize = stepsize, prestepsize2 = 0, f2pre = 0;
		if (maxpref - f2 < -LS_alpha * initialslope * stepsize)
		{
			stepsize = -initialslope * prestepsize * prestepsize / 2 / (f2 - f1 - initialslope * prestepsize);
			stepsize = (stepsize < LS_ratio2 * prestepsize) ? stepsize : LS_ratio2 * prestepsize;
			stepsize = (stepsize > LS_ratio1 * prestepsize) ? stepsize : LS_ratio1 * prestepsize;
			f2pre = f2;
			prestepsize2 = prestepsize;
			f2 = h();
			prestepsize = stepsize;
		}

		realdp a11, a12, a21, a22, b1, b2, c, a = 0, b = 0;
		while (maxpref - f2 < -LS_alpha * initialslope * stepsize)
		{
			a11 = static_cast<realdp> (1) / prestepsize / prestepsize;
			a12 = static_cast<realdp> (- 1) / prestepsize2 / prestepsize2;
			a21 = -prestepsize2 / prestepsize / prestepsize;
			a22 = prestepsize / prestepsize2 / prestepsize2;
			b1 = f2 - f1 - initialslope * prestepsize;
			b2 = f2pre - f1 - initialslope * prestepsize2;
			c = prestepsize - prestepsize2;
			a = (a11 * b1 + a12 * b2) / c;
			b = (a21 * b1 + a22 * b2) / c;
			stepsize = (-b + sqrt(b * b - 3 * a * initialslope)) / 3 / a;
			stepsize = (stepsize < LS_ratio2 * prestepsize) ? stepsize : LS_ratio2 * prestepsize;
			stepsize = (stepsize > LS_ratio1 * prestepsize) ? stepsize : LS_ratio1 * prestepsize;
			if (stepsize < Minstepsize)
			{
				if (Verbose > FINALRESULT)
				{
					printf("Warning: step size reaches the minimum: %3.2e!\n", Minstepsize);
				}
				LSstatus = LSSM_MINSTEPSIZE;
				break;
			}
			f2pre = f2;
			prestepsize2 = prestepsize;
			f2 = h();
			prestepsize = stepsize;
		}
        Prob->Grad(x2, &gf2); ng++;
        newslope = INFINITY;
	};

	void SolversSMLS::LinesearchExact(void)
	{
		realdp a1, B, t, dir, s, y, sgf1, sgf2, sf1, sf2 = 0, LS_ratio = LS_ratio1;
        realdp minstep = 100 * std::sqrt(std::numeric_limits<realdp>::epsilon());
		LSstatus = LSSM_SUCCESS;
		t = 1;                /* initial step size */
		a1 = 0;				  /* initial iterate is 0 */
		sgf1 = initialslope;  /* gradient at initial iterate */
		B = fabs(sgf1 / initiallength);   /* initial Hessian approximation */
		sf1 = f1;
        sgf2 = sgf1;
        integer times = 0;
		while (times < 20)
		{   /* quasi-Newton for optimizing the scalar function */
			dir = - sgf1 / B; /* obtain direction - B^{-1} grad */
			/* start Armijo-Goldstein line search */
			t = 1;         /* initial step size is 1 */
			stepsize = a1 + t * dir; sf2 = h();   /* evaluete the scalar function at (a1 + t * dir) */
			while (sf2 > sf1 + LS_alpha * t * sgf1 * dir)
			{
				t *= LS_ratio;
				stepsize = a1 + t * dir; sf2 = h();
				if (t * LS_ratio < minstep)
				{
					sgf2 = dh();
					if (fabs(sgf2 / initialslope) > 1e-1 && sf2 > f1 + LS_alpha * stepsize * ngf1)
						LSstatus = LSSM_NONEXACT;
					break;
				}
			}

			if (LSstatus == LSSM_NONEXACT)
				break;
            
			s = t * dir;
			if (fabs(s) < minstep)
			{
				sgf2 = dh();
				break;
			}
			sgf2 = dh();
            
#ifdef SINGLE_PRECISION
			if (fabs(sgf2 / initialslope) < 1e-2)
				break;
#else
            if (fabs(sgf2 / initialslope) < 1e-6)
                break;
#endif
			y = sgf2 - sgf1;
			B = (y / s > 0) ? y / s : B;     /* for scalar function, Hessian approximation is uniquely defined by s and y. */
			/* update */
			a1 += t * dir;
			sgf1 = sgf2;
			sf1 = sf2;
            times++;
		}

		if (stepsize <= Minstepsize)
		{
			LSstatus = LSSM_MINSTEPSIZE;
		}
		if (stepsize >= Maxstepsize)
		{
			LSstatus = LSSM_MAXSTEPSIZE;
		}
		f2 = sf2;
		newslope = sgf2;
	};

	realdp SolversSMLS::h(void)
	{
        Mani->ScalarTimesVector(x1, stepsize, eta1, &eta2);
        Mani->Retraction(x1, eta2, &x2); nR++;
        nf++;
        return Prob->f(x2);
	};

	realdp SolversSMLS::dh(void)
	{
		Prob->Grad(x2, &gf2); ng++;
		nV++;
        Vector diffeta1(Mani->GetEMPTY());
        Mani->DiffRetraction(x1, eta2, x2, eta1, &diffeta1, true);
        realdp tmp = Mani->Metric(x2, gf2, diffeta1);
		return tmp;
	};

	/* Algorithm 3.5 in NW06 */
	void SolversSMLS::LinesearchStrongWolfe(void)
	{
		realdp prestepsize = 0, fpre = std::numeric_limits<realdp>::max(), newslopepre = initialslope;
        realdp tstepsize = 0;
		LSstatus = LSSM_SUCCESS;
		while (1)
		{
			f2 = h();
			if (f2 > f1 + LS_alpha * stepsize * initialslope || f2 >= fpre)
			{
				Zoom(prestepsize, fpre, newslopepre, stepsize, f2);
				return;
			}
			newslope = dh();
			if (fabs(newslope) <= -LS_beta * initialslope)
			{
				return;
			}
			if (newslope >= 0)
			{
				Zoom(stepsize, f2, newslope, prestepsize, fpre);
				return;
			}
			prestepsize = stepsize;
			fpre = f2;
			newslopepre = newslope;
			if (stepsize == Maxstepsize)
			{
				LSstatus = LSSM_MAXSTEPSIZE;
				return;
			}
            tstepsize = - initialslope * stepsize * stepsize / 2 / (f2 - f1 - initialslope * stepsize);
            stepsize = (1.1 * stepsize < tstepsize) ? tstepsize : 1.1 * stepsize;
			stepsize = (stepsize < Maxstepsize) ? stepsize : Maxstepsize;
		}
	};

	void SolversSMLS::Zoom(realdp x1, realdp fx1, realdp slopex1, realdp x2, realdp fx2)
	{
		realdp xdiff, xincr, xlo = x1, xhi = x2, fxlo = fx1, fxhi = fx2, xloslope = slopex1;
		integer times = 0;
		while (1)
		{
			xdiff = (xhi - xlo);
			xincr = -xloslope * xdiff * xdiff / 2 / (fxhi - (fxlo + xloslope * xdiff));
            xincr = (xincr < xdiff * LS_ratio1) ? xdiff * LS_ratio1 : xincr;
            xincr = (xincr < xdiff * LS_ratio2) ? xincr : xdiff * LS_ratio2;
			stepsize = xlo + xincr;
			times++;
			if (times >= 10)
			{
				LSstatus = LSSM_LSERROR;
				return;
			}
			f2 = h();
			if (f2 > f1 + LS_alpha * stepsize * initialslope || f2 >= fxlo)
			{
				xhi = stepsize;
				fxhi = f2;
			}
			else
			{
				newslope = dh();
				if (fabs(newslope) <= -LS_beta * initialslope)
				{
					return;
				}
				if (newslope * (xhi - xlo) >= 0)
				{
					xhi = xlo;
					fxhi = fxlo;
				}
				xlo = stepsize;
				fxlo = f2;
				xloslope = newslope;
			}
			if (stepsize <= Minstepsize)
			{
				LSstatus = LSSM_MINSTEPSIZE;
				return;
			}
		}
	};

	void SolversSMLS::CheckParams(void)
	{
		SolversSM::CheckParams();

		std::string LSALGOnames[LSALGOSMLENGTH] = { "LSSM_ARMIJO", "LSSM_WOLFE", "LSSM_STRONGWOLFE", "LSSM_EXACT", "LSSM_INPUTFUN" };
		std::string INITSTEPnames[INITSTEPSIZESETSMLENGTH] = { "LSSM_ONESTEP", "LSSM_BBSTEP", "LSSM_QUADINT", "LSSM_QUADINTMOD", "LSSM_EXTRBBSTEP" };

		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		printf("LINE SEARCH TYPE METHODS PARAMETERS:\n");
		status = (LineSearch_LS >= 0 && LineSearch_LS < LSALGOSMLENGTH) ? YES : NO;
		printf("LineSearch_LS :%15s[%s],\t", LSALGOnames[LineSearch_LS].c_str(), status);
        status = (InitSteptype >= 0 && InitSteptype < INITSTEPSIZESETSMLENGTH) ? YES : NO;
        printf("InitSteptype  :%15s[%s],\n", INITSTEPnames[InitSteptype].c_str(), status);
		status = (LS_alpha > 0 && LS_alpha < 0.5) ? YES : NO;
		printf("LS_alpha      :%15g[%s],\t", LS_alpha, status);
        status = (LS_beta > LS_alpha && LS_beta < 1) ? YES : NO;
        printf("LS_beta       :%15g[%s],\n", LS_beta, status);
        status = (LS_ratio1 > 0 && LS_ratio1 <= LS_ratio2) ? YES : NO;
        printf("LS_ratio1     :%15g[%s],\t", LS_ratio1, status);
        status = (LS_ratio2 >= LS_ratio1 && LS_ratio2 < 1) ? YES : NO;
        printf("LS_ratio2     :%15g[%s],\n", LS_ratio2, status);
		status = (Initstepsize > 0) ? YES : NO;
		printf("Initstepsize  :%15g[%s],\t", Initstepsize, status);
        status = YES;
        printf("Finalstepsize :%15g[%s],\n", Finalstepsize, status);
		status = (Minstepsize > 0 && Minstepsize <= Maxstepsize) ? YES : NO;
		printf("Minstepsize   :%15g[%s],\t", Minstepsize, status);
		status = (Maxstepsize > 0 && Maxstepsize >= Minstepsize) ? YES : NO;
		printf("Maxstepsize   :%15g[%s],\n", Maxstepsize, status);
		status = (Num_pre_funs >= 1) ? YES : NO;
		printf("Num_pre_funs  :%15d[%s],\t", Num_pre_funs, status);
		status = (Num_pre_BB >= 0) ? YES : NO;
		printf("Num_pre_BB    :%15d[%s],\n", Num_pre_BB, status);
        status = ((BBratio >= 0 && BBratio <= 1) || BBratio == -1) ? YES : NO;
        printf("BBratio       :%15g[%s],\t", BBratio, status);
		status = (PreFunsAccuracy >= 0) ? YES : NO;
		printf("PreFunsAccuracy:%14g[%s],\n", PreFunsAccuracy, status);
		if (LineSearch_LS == LSSM_INPUTFUN)
		{
			status = YES;
			printf("IsPureLSInput :%15d[%s]\n", IsPureLSInput, status);
		}
	};

	void SolversSMLS::UpdateData(void)
	{
	};

    void SolversSMLS::PrintInfo(void)
    {
        printf("i:%d,f:%.3e,df/f:%.3e,", iter, f2, ((f1 - f2) / std::fabs(f2)));

        printf("|gf|:%.3e,t0:%.2e,t:%.2e,s0:%.2e,s:%.2e,time:%.2g,", ngf2, initiallength, stepsize, initialslope, newslope, static_cast<realdp>(getTickCount() - starttime) / CLK_PS);

        printf("nf:%d,ng:%d,", nf, ng);
        
        if (nH != 0)
            printf("nH:%d,", nH);
        
        printf("nR:%d,", nR);
        
        if (nV != 0)
            printf("nV(nVp):%d(%d),", nV, nVp);
        
        printf("\n");
    };

    void SolversSMLS::PrintFinalInfo(void)
    {
        printf("i:%d,f:%.3e,", iter, f2);

        printf("|gf|:%.3e,|gf|/|gf0|:%.3e,time:%.2g,", ngf2, ngf2/ngf0, static_cast<realdp>(getTickCount() - starttime) / CLK_PS);
        
        printf("nf:%d,ng:%d,", nf, ng);
        
        if (nH != 0)
            printf("nH:%d,", nH);
        
        printf("nR:%d,", nR);
        
        if (nV != 0)
            printf("nV(nVp):%d(%d),", nV, nVp);
        
        printf("\n");
    };

	void SolversSMLS::SetDefaultParams()
	{
		SolversSM::SetDefaultParams();
		LineSearch_LS = LSSM_ARMIJO;
		LinesearchInput = nullptr;
		Num_pre_funs = 1;
		LS_alpha = static_cast<realdp> (1e-4);
		LS_beta = static_cast<realdp> (0.999);
		Minstepsize = std::numeric_limits<realdp>::epsilon();
		Maxstepsize = static_cast<realdp> (1e10);
		LS_ratio1 = static_cast<realdp> (0.1);
		LS_ratio2 = static_cast<realdp> (0.9);
		Initstepsize = 1;
		Finalstepsize = -1;
		Num_pre_BB = 0;
		BBratio = static_cast<realdp> (1);
		PreFunsAccuracy = 1e6;
		IsPureLSInput = false;
		LSstatusSetnames = new std::string[LSSTATUSSETSMLENGTH];
		LSstatusSetnames[LSSM_NOCURVATURE].assign("LSSM_NOCURVATURE");
		LSstatusSetnames[LSSM_MINSTEPSIZE].assign("LSSM_MINSTEPSIZE");
		LSstatusSetnames[LSSM_MAXSTEPSIZE].assign("LSSM_MAXSTEPSIZE");
		LSstatusSetnames[LSSM_NONEXACT].assign("LSSM_NONEXACT");
		LSstatusSetnames[LSSM_LSERROR].assign("LSSM_LSERROR");
		LSstatusSetnames[LSSM_SUCCESS].assign("LSSM_SUCCESS");
	};

	SolversSMLS::~SolversSMLS(void)
	{
		delete[] LSstatusSetnames;
	};

	/* Algorithm 6.3.1MOD in [DS83] */
	void SolversSMLS::LinesearchWolfe(void)
	{
		realdp prestepsize = 0, f2pre = 0;
		realdp stepsizelo, stepsizediff, stepsizeincr, steptemp;
		realdp flo, fhi;
		integer times = 0;
		realdp n1, n2, n3, n4, a, b, disc;
		LSstatus = LSSM_SUCCESS;
		while (1)
		{
			f2 = h();
			if (f2 <= f1 + LS_alpha * stepsize * initialslope)
			{
				newslope = dh();
				if (newslope < LS_beta * initialslope)
				{
					times = 0;
					while (f2 <= f1 + LS_alpha * stepsize * initialslope && newslope < LS_beta * initialslope && stepsize < Maxstepsize)
					{
						prestepsize = stepsize;
						f2pre = f2;
						stepsize = (2 * stepsize < Maxstepsize) ? 2 * stepsize : Maxstepsize;
						f2 = h();
						if (f2 <= f1 + LS_alpha * stepsize * initialslope)
							newslope = dh();
						times++;
						if (times > 10)
							break;
					}
					if (stepsize >= Maxstepsize)
					{
						Prob->Grad(x2, &gf2); ng++;
						LSstatus = LSSM_MAXSTEPSIZE;
						return;
					}
					if (stepsize != initiallength && f2 > f1 + LS_alpha * stepsize * initialslope)
					{
						stepsizelo = (stepsize < prestepsize) ? stepsize : prestepsize;
						stepsizediff = fabs(prestepsize - stepsize);
						if (stepsize < prestepsize)
						{
							flo = f2;
							fhi = f2pre;
						}
						else
						{
							flo = f2pre;
							fhi = f2;
						}
						times = 0;
						while ((f2 > f1 + LS_alpha * stepsize * initialslope || newslope < LS_beta * initialslope) && stepsizediff >= Minstepsize)
						{
							stepsizeincr = -newslope * stepsizediff * stepsizediff / 2 / (fhi - (flo + newslope * stepsizediff));
							stepsizeincr = (stepsizeincr < static_cast<realdp> (0.2) * stepsizediff) ? static_cast<realdp> (0.2) * stepsizediff : stepsizeincr;
							stepsize = stepsizelo + stepsizeincr;
							f2 = h();
							if (f2 > f1 + LS_alpha * stepsize * initialslope)
							{
								stepsizediff = stepsizeincr;
								fhi = f2;
							}
							else
							{
								newslope = dh();
								if (newslope < LS_beta * initialslope)
								{
									stepsizelo = stepsize;
									stepsizediff -= stepsizeincr;
									flo = f2;
								}
							}
							times++;
							if (times > 10)
								break;
						}
						if (newslope < LS_beta * initialslope)
						{
							f2 = h();
							newslope = dh();
							LSstatus = LSSM_NOCURVATURE;
							return;
						}
					}
				}
				LSstatus = LSSM_SUCCESS;
				return;
			}
			if (stepsize <= Minstepsize)
			{
				stepsize = Minstepsize;
				f2 = h();
				newslope = dh();
				LSstatus = LSSM_MINSTEPSIZE;
				return;
			}
			else
			{
				if (stepsize == initiallength)
					steptemp = -initialslope * initiallength * initiallength / 2 / (f2 - f1 - initialslope * initiallength);
				else
				{
					n1 = static_cast<realdp> (1) / stepsize / stepsize;
					n2 = static_cast<realdp> (1) / prestepsize / prestepsize;
					n3 = (f2 - f1 - stepsize * initialslope) / (stepsize - prestepsize);
					n4 = (f2pre - f1 - prestepsize * initialslope) / (stepsize - prestepsize);
					a = n1 * n3 - n2 * n4;
					b = -prestepsize * n1 * n3 + stepsize * n2 * n4;
					disc = b * b - 3 * a * initialslope;
					disc = (disc < 0) ? 0 : disc;
					if (fabs(a) < 1e-10)
						steptemp = -initialslope / 2 / b;
					else
						steptemp = (-b + sqrt(disc)) / 3 / a;
					steptemp = (steptemp > 0.5 * stepsize) ? static_cast<realdp> (0.5) * stepsize : steptemp;
				}
				prestepsize = stepsize;
				f2pre = f2;
				stepsize = (steptemp <= static_cast<realdp> (1e-2) * stepsize) ? static_cast<realdp> (1e-2) * stepsize : steptemp;
			}
		}
	};

	void SolversSMLS::SetParams(PARAMSMAP params)
	{
		SolversSM::SetParams(params);
		PARAMSMAP::iterator iter;
		for (iter = params.begin(); iter != params.end(); iter++)
		{
			if (iter->first == static_cast<std::string> ("LineSearch_LS"))
			{
				LineSearch_LS = static_cast<LSAlgoSM> (static_cast<integer> (iter->second));
			}
			else
			if (iter->first == static_cast<std::string> ("LS_alpha"))
			{
				LS_alpha = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("LS_beta"))
			{
				LS_beta = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("Minstepsize"))
			{
				Minstepsize = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("Maxstepsize"))
			{
				Maxstepsize = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("LS_ratio1"))
			{
				LS_ratio1 = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("LS_ratio2"))
			{
				LS_ratio2 = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("Initstepsize"))
			{
				Initstepsize = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("Finalstepsize"))
			{
				Finalstepsize = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("Num_pre_funs"))
			{
				Num_pre_funs = static_cast<integer> (iter->second);
			}
			else
			if (iter->first == static_cast<std::string> ("InitSteptype"))
			{
				InitSteptype = static_cast<InitStepsizeSetSM> (static_cast<integer> (iter->second));
			}
			else
			if (iter->first == static_cast<std::string> ("IsPureLSInput"))
			{
				IsPureLSInput = ((static_cast<integer> (iter->second)) != 0);
			}
			else
			if (iter->first == static_cast<std::string> ("PreFunsAccuracy"))
			{
				PreFunsAccuracy = iter->second;
			}
            else
            if (iter->first == static_cast<std::string> ("Num_pre_BB"))
            {
                Num_pre_BB = static_cast<integer> (iter->second);
            }
            else
            if (iter->first == static_cast<std::string> ("BBratio"))
            {
                BBratio = static_cast<realdp> (iter->second);
            }
		}
	};
}; /*end of ROPTLIB namespace*/
