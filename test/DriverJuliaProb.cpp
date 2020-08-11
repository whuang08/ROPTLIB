
#include "test/DriverJuliaProb.h"

#ifdef DRIVERJULIAPROB

using namespace ROPTLIB;

//double *DriverJuliaProb(const char *fname, const char *gfname, const char *hfname, const char *isstopped, const char *LSinput, const char *solvername, double *paramsvalues, long int lengthSParams)
double *DriverJuliaProb(const char *fname, const char *gfname, const char *hfname, const char *isstopped, const char *LSinput, /*Handles*/
                        const char *solvername, double *paramsvalues, long long lengthSParams, /*Sparams*/
                        const char *maninames, long long numoftypes, long long *numofmani, long long *paramset, long long *ms, long long *ns, long long *ps, long long IsCheckParams,
                        long long inHasHHR, double *X0, integer length_X0)
{
    FunHandles *Handles = new FunHandles();
    SolverParams *Sparams = new SolverParams();
    ManiParams *Mparams = new ManiParams();

    Handles->fname = fname;
    Handles->gfname = gfname;
    Handles->hfname = hfname;
    Handles->isstopped = isstopped;
    Handles->LinesearchInput = LSinput;
    
    Sparams->name = solvername;
    Sparams->IsCheckParams = paramsvalues[0];
    Sparams->IsCheckGradHess = paramsvalues[1];
    Sparams->Stop_Criterion = paramsvalues[2];
    Sparams->Tolerance = paramsvalues[3];
    Sparams->Diffx = paramsvalues[4];
    Sparams->NumExtraGF = paramsvalues[5];
    Sparams->TimeBound = paramsvalues[6];
    Sparams->Min_Iteration = paramsvalues[7];
    Sparams->Max_Iteration = paramsvalues[8];
    Sparams->OutputGap = paramsvalues[9];
    Sparams->Verbose = paramsvalues[10];
    Sparams->isconvex = paramsvalues[11];
    Sparams->nu = paramsvalues[12];
    Sparams->mu = paramsvalues[13];
    Sparams->LengthSY = paramsvalues[14];
    Sparams->lambdaLower = paramsvalues[15];
    Sparams->lambdaUpper = paramsvalues[16];
    Sparams->LineSearch_LS = paramsvalues[17];
    Sparams->IsPureLSInput = paramsvalues[18];
    Sparams->LS_alpha = paramsvalues[19];
    Sparams->LS_beta = paramsvalues[20];
    Sparams->Minstepsize = paramsvalues[21];
    Sparams->Maxstepsize = paramsvalues[22];
    Sparams->LS_ratio1 = paramsvalues[23];
    Sparams->LS_ratio2 = paramsvalues[24];
    Sparams->Initstepsize = paramsvalues[25];
    Sparams->Accuracy = paramsvalues[26];
    Sparams->Finalstepsize = paramsvalues[27];
    Sparams->Num_pre_funs = paramsvalues[28];
    Sparams->InitSteptype = paramsvalues[29];
    Sparams->Acceptence_Rho = paramsvalues[30];
    Sparams->Shrinked_tau = paramsvalues[31];
    Sparams->Magnified_tau = paramsvalues[32];
    Sparams->minimum_Delta = paramsvalues[33];
    Sparams->maximum_Delta = paramsvalues[34];
    Sparams->useRand = paramsvalues[35];
    Sparams->Max_Inner_Iter = paramsvalues[36];
    Sparams->Min_Inner_Iter = paramsvalues[37];
    Sparams->theta = paramsvalues[38];
    Sparams->kappa = paramsvalues[39];
    Sparams->initial_Delta = paramsvalues[40];
    Sparams->Eps = paramsvalues[41];
    Sparams->Theta_eps = paramsvalues[42];
    Sparams->Min_Eps = paramsvalues[43];
    Sparams->Del = paramsvalues[44];
    Sparams->Theta_del = paramsvalues[45];
    const char **Maninames = new const char*[numoftypes];
    char *Maninamestr = new char[1000];
    strcpy(Maninamestr, maninames);
    integer idx = 0;
    Maninames[idx] = strtok(Maninamestr, " ,.-");
//    printf("hhh\n");//---
//    std::cout << "numoftypes:" << numoftypes << std::endl;//---
    std::cout << Maninames[idx] << std::endl;//---
    while(Maninames[idx] != nullptr)
    {
//        std::cout << "idx:" << idx << std::endl;//---
        idx++;
        if(idx == numoftypes)
            break;
        Maninames[idx] = strtok(nullptr, " ,.-");
        std::cout << Maninames[idx] << std::endl;//---
    }
    Mparams->name = Maninames;
    Mparams->m = ms;
    Mparams->n = ns;
    Mparams->p = ps;
    Mparams->numoftypes = numoftypes;
    Mparams->numofmani = numofmani;
    Mparams->paramset = paramset;
    Mparams->IsCheckParams = IsCheckParams;
    
//    std::cout << Handles->fname << std::endl;//---
//    std::cout << Handles->gfname << std::endl;//---
//    std::cout << Handles->hfname << std::endl;//---
//    std::cout << Handles->isstopped << std::endl;//---
//    std::cout << Handles->LinesearchInput << std::endl;//---
//    std::cout << Sparams->name << std::endl;//---
//    std::cout << maninames << std::endl;//---
    
    // Initialization of Julia is not necessary.
    // If this code is run in C++ environment,  then calling julia requires initialization of Julia.
    // However, this C++ code is run in Julia. This implies Julia has been run when this code is called.
    // Therefore, it is not necessary to initialize Julia.
    //	jl_init(JULIA_DIR)

    // Obtain manifold and iterate structure
    Manifold *domain = nullptr, **manifolds = nullptr;
    Variable initialX;
    integer *powsinterval, numoftype, numoftotal;

    if (!ParseManiParams(Mparams, manifolds, numoftype, powsinterval))
    {
        std::cout << "Parsing ManiParams fails." << std::endl;
        exit(EXIT_FAILURE);
    }

    numoftotal = powsinterval[numoftype];
    if(numoftotal > 1)
        domain = new MultiManifolds(manifolds, numoftype, powsinterval);
    else
        domain = manifolds[0];
    
    domain->SetHasHHR(inHasHHR != 0);
    if(Mparams->IsCheckParams != 0)
        domain->CheckParams();

    initialX = domain->RandominManifold();
//    std::cout << "h1" << std::endl;//---
    // initialize the initial iterate
    if (length_X0 != 0)
    {
        double *initXptr = initialX.ObtainWriteEntireData();
        if(length_X0 != initialX.Getlength())
        {
            std::cout << "Error: The initial iterate does not have correct size: " << length_X0 << "!=" << initialX.Getlength() << "!" << std::endl;
            exit(EXIT_FAILURE);
        }
        dcopy_(&length_X0, X0, &GLOBAL::IONE, initXptr, &GLOBAL::IONE);
    }
    
//    std::cout << "h2" << std::endl;//---
    // Define the problem
    jl_function_t *func = jl_get_function(jl_main_module, Handles->fname);
    jl_function_t *gfunc = jl_get_function(jl_main_module, Handles->gfname);
    jl_function_t *hfunc = jl_get_function(jl_main_module, Handles->hfname);
    
//    std::cout << (long) func << std::endl;//---
//    std::cout << (long) gfunc << std::endl;//---
//    std::cout << (long) hfunc << std::endl;//---
    
//    std::cout << "h3" << std::endl;//---
    Problem *Prob = new juliaProblem(func, gfunc, hfunc);
//    std::cout << "h31" << std::endl;//---
    Prob->SetDomain(domain);
    
//    std::cout << "h4" << std::endl;//---
    double *result = ParseSolverParamsAndOptimizing(Sparams, Handles, Prob, &initialX);
//    std::cout << "h5" << std::endl;//---
    delete Prob;
//    std::cout << "h6" << std::endl;//---
    if(numoftotal > 1)
        delete domain;
    
//    std::cout << "h7" << std::endl;//---
    for (integer i = 0; i < numoftype; i++)
    {
        delete manifolds[i];
    }
//    std::cout << "h8" << std::endl;//---
    delete[] manifolds;
//    std::cout << "h9" << std::endl;//---
    delete[] powsinterval;
    
//    std::cout << "h10" << std::endl;//---
    delete Handles;
//    std::cout << "h11" << std::endl;//---
    delete Sparams;
//    std::cout << "h12" << std::endl;//---
    delete Mparams;
//    std::cout << "h13" << std::endl;//---
    delete [] Maninames;
//    std::cout << "h14" << std::endl;//---
    delete [] Maninamestr;
//    std::cout << "h15" << std::endl;//---
    
    return result;

//    It is not necessary to close Julia. The reason is the same as that in "jl_init"
//    jl_atexit_hook(0);
//    return;
};

namespace RJULIA{
    jl_function_t *isstopped = nullptr;
    /*This function defines the stopping criterion that may be used in the C++ solver*/
    bool juliaInnerStop(const Variable &x, const Vector &funSeries, integer lengthSeries, realdp ngf, realdp ngf0, const Problem *prob, const Solvers *solver)
    {
        jl_value_t *args[4] = {nullptr};
        jl_value_t* array_type = jl_apply_array_type((jl_value_t *) jl_float64_type, 1);
        const double *xptr = x.ObtainReadData();
        jl_array_t *arrx = jl_ptr_to_array_1d(array_type, const_cast<double *> (xptr), x.Getlength(), 0);
        const double *funptr = funSeries.ObtainReadData();
        jl_array_t *arrfuns = jl_ptr_to_array_1d(array_type, const_cast<double *> (funptr), lengthSeries, 0);
        args[0] = (jl_value_t *) arrx;
        args[1] = (jl_value_t *) arrfuns;
        args[2] = jl_box_float64(ngf);
        args[3] = jl_box_float64(ngf0);

        jl_value_t *retx = jl_call(isstopped, args, 4);

        if(jl_is_int64(retx))
        {
            integer result = jl_unbox_int64(retx);
            return (result != 0);
        }
        if(jl_is_bool(retx))
        {
            return jl_unbox_bool(retx);
        }

        std::cout << "Error: Function isstopped must return an integer or a bool!" << std::endl;
        exit(EXIT_FAILURE);
    };

    jl_function_t *LinesearchInput = nullptr;
    /*This function defines the line search algorithm that may be used in the C++ solver*/
    double juliaLinesearchInput(integer iter, const Variable &x1, const Vector &eta1, realdp initialstepsize, realdp initialslope, const Problem *prob, const Solvers *solver)
    {
        jl_value_t *args[5] = {nullptr};
        jl_value_t* array_type = jl_apply_array_type((jl_value_t *) jl_float64_type, 1);
        const double *xptr = x1.ObtainReadData();
        jl_array_t *arrx = jl_ptr_to_array_1d(array_type, const_cast<double *> (xptr), x1.Getlength(), 0);
        const double *etaptr = eta1.ObtainReadData();
        jl_array_t *arreta = jl_ptr_to_array_1d(array_type, const_cast<double *> (etaptr), eta1.Getlength(), 0);
        args[0] = (jl_value_t *) arrx;
        args[1] = (jl_value_t *) arreta;
        args[2] = jl_box_float64(initialstepsize);
		args[3] = jl_box_float64(initialslope);
		args[4] = jl_box_int64(iter);

        jl_value_t *retx = jl_call(LinesearchInput, args, 5);

        double result = jl_unbox_float64(retx);
        return result;
    }
};

double *ParseSolverParamsAndOptimizing(struct SolverParams *Sparams, struct FunHandles *Handles, Problem *Prob, Variable *initialX)
{
    PARAMSMAP params;

//    Add parameters
    if(Sparams->Stop_Criterion != -1)
        params.insert(std::pair<std::string, double>("Stop_Criterion", Sparams->Stop_Criterion));
    if(Sparams->Tolerance != -1)
        params.insert(std::pair<std::string, double>("Tolerance", Sparams->Tolerance));
    if(Sparams->Diffx != -1)
        params.insert(std::pair<std::string, double>("Diffx", Sparams->Diffx));
    if(Sparams->NumExtraGF != -1)
        params.insert(std::pair<std::string, double>("NumExtraGF", Sparams->NumExtraGF));
    if(Sparams->TimeBound != -1)
        params.insert(std::pair<std::string, double>("TimeBound", Sparams->TimeBound));
    if(Sparams->Min_Iteration != -1)
        params.insert(std::pair<std::string, double>("Min_Iteration", Sparams->Min_Iteration));
    if(Sparams->Max_Iteration != -1)
        params.insert(std::pair<std::string, double>("Max_Iteration", Sparams->Max_Iteration));
    if(Sparams->OutputGap != -1)
        params.insert(std::pair<std::string, double>("OutputGap", Sparams->OutputGap));
    if(Sparams->Verbose != -1)
        params.insert(std::pair<std::string, double>("Verbose", Sparams->Verbose));
    if(Sparams->isconvex != -1)
        params.insert(std::pair<std::string, double>("isconvex", Sparams->isconvex));
    if(Sparams->nu != -1)
        params.insert(std::pair<std::string, double>("nu", Sparams->nu));
    if(Sparams->mu != -1)
        params.insert(std::pair<std::string, double>("mu", Sparams->mu));
    if(Sparams->LengthSY != -1)
        params.insert(std::pair<std::string, double>("LengthSY", Sparams->LengthSY));
    if(Sparams->lambdaLower != -1)
        params.insert(std::pair<std::string, double>("lambdaLower", Sparams->lambdaLower));
    if(Sparams->lambdaUpper != -1)
        params.insert(std::pair<std::string, double>("lambdaUpper", Sparams->lambdaUpper));
    if(Sparams->LineSearch_LS != -1)
        params.insert(std::pair<std::string, double>("LineSearch_LS", Sparams->LineSearch_LS));
    if(Sparams->IsPureLSInput != -1)
        params.insert(std::pair<std::string, double>("IsPureLSInput", Sparams->IsPureLSInput));
    if(Sparams->LS_alpha != -1)
        params.insert(std::pair<std::string, double>("LS_alpha", Sparams->LS_alpha));
    if(Sparams->LS_beta != -1)
        params.insert(std::pair<std::string, double>("LS_beta", Sparams->LS_beta));
    if(Sparams->Minstepsize != -1)
        params.insert(std::pair<std::string, double>("Minstepsize", Sparams->Minstepsize));
    if(Sparams->Maxstepsize != -1)
        params.insert(std::pair<std::string, double>("Maxstepsize", Sparams->Maxstepsize));
    if(Sparams->LS_ratio1 != -1)
        params.insert(std::pair<std::string, double>("LS_ratio1", Sparams->LS_ratio1));
    if(Sparams->LS_ratio2 != -1)
        params.insert(std::pair<std::string, double>("LS_ratio2", Sparams->LS_ratio2));
    if(Sparams->Initstepsize != -1)
        params.insert(std::pair<std::string, double>("Initstepsize", Sparams->Initstepsize));
    if(Sparams->Accuracy != -1)
        params.insert(std::pair<std::string, double>("Accuracy", Sparams->Accuracy));
    if(Sparams->Finalstepsize != -1)
        params.insert(std::pair<std::string, double>("Finalstepsize", Sparams->Finalstepsize));
    if(Sparams->Num_pre_funs != -1)
        params.insert(std::pair<std::string, double>("Num_pre_funs", Sparams->Num_pre_funs));
    if(Sparams->InitSteptype != -1)
        params.insert(std::pair<std::string, double>("InitSteptype", Sparams->InitSteptype));
    if(Sparams->Acceptence_Rho != -1)
        params.insert(std::pair<std::string, double>("Acceptence_Rho", Sparams->Acceptence_Rho));
    if(Sparams->Shrinked_tau != -1)
        params.insert(std::pair<std::string, double>("Shrinked_tau", Sparams->Shrinked_tau));
    if(Sparams->Magnified_tau != -1)
        params.insert(std::pair<std::string, double>("Magnified_tau", Sparams->Magnified_tau));
    if(Sparams->minimum_Delta != -1)
        params.insert(std::pair<std::string, double>("minimum_Delta", Sparams->minimum_Delta));
    if(Sparams->maximum_Delta != -1)
        params.insert(std::pair<std::string, double>("maximum_Delta", Sparams->maximum_Delta));
    if(Sparams->useRand != -1)
        params.insert(std::pair<std::string, double>("useRand", Sparams->useRand));
    if(Sparams->Max_Inner_Iter != -1)
        params.insert(std::pair<std::string, double>("Max_Inner_Iter", Sparams->Max_Inner_Iter));
    if(Sparams->Min_Inner_Iter != -1)
        params.insert(std::pair<std::string, double>("Min_Inner_Iter", Sparams->Min_Inner_Iter));
    if(Sparams->theta != -1)
        params.insert(std::pair<std::string, double>("theta", Sparams->theta));
    if(Sparams->kappa != -1)
        params.insert(std::pair<std::string, double>("kappa", Sparams->kappa));
    if(Sparams->initial_Delta != -1)
        params.insert(std::pair<std::string, double>("initial_Delta", Sparams->initial_Delta));
    if(Sparams->Eps != -1)
        params.insert(std::pair<std::string, double>("Eps", Sparams->Eps));
    if(Sparams->Theta_eps != -1)
        params.insert(std::pair<std::string, double>("Theta_eps", Sparams->Theta_eps));
    if(Sparams->Min_Eps != -1)
        params.insert(std::pair<std::string, double>("Min_Eps", Sparams->Min_Eps));
    if(Sparams->Del != -1)
        params.insert(std::pair<std::string, double>("Del", Sparams->Del));
    if(Sparams->Theta_del != -1)
        params.insert(std::pair<std::string, double>("Theta_del", Sparams->Theta_del));

    std::cout << "h41" << std::endl;//---
    std::string stdmethodname = Sparams->name;
    Solvers *solver;
    if (stdmethodname == "RSD")
    {
        solver = new RSD(Prob, initialX);
    }
    else
    if (stdmethodname == "RNewton")
    {
		solver = new RNewton(Prob, initialX);
    }
    else
    if (stdmethodname == "RCG")
    {
		solver = new RCG(Prob, initialX);
    }
    else
    if (stdmethodname == "RBroydenFamily")
    {
		solver = new RBroydenFamily(Prob, initialX, nullptr);
    }
    else
    if (stdmethodname == "RWRBFGS")
    {
		solver = new RWRBFGS(Prob, initialX, nullptr);
    }
    else
    if (stdmethodname == "RBFGS")
    {
		solver = new RBFGS(Prob, initialX, nullptr);
    }
    else
    if (stdmethodname == "RBFGSSub")
    {
		solver = new RBFGSSub(Prob, initialX, nullptr);
    }
    else
    if (stdmethodname == "LRBFGSSub")
    {
		solver = new LRBFGSSub(Prob, initialX);
    }
    else
    if (stdmethodname == "RGS")
    {
		solver = new RGS(Prob, initialX);
    }
    else
    if (stdmethodname == "LRBFGS")
    {
		solver = new LRBFGS(Prob, initialX);
    }
    else
    if (stdmethodname == "RTRSD")
    {
		solver = new RTRSD(Prob, initialX);
    }
    else
    if (stdmethodname == "RTRNewton")
    {
		solver = new RTRNewton(Prob, initialX);
    }
    else
    if (stdmethodname == "RTRSR1")
    {
		solver = new RTRSR1(Prob, initialX, nullptr);
    }
    else
    if (stdmethodname == "LRTRSR1")
    {
		solver = new LRTRSR1(Prob, initialX);
    }
    else
    if (stdmethodname == "ManPG")
    {
        solver = new ManPG(Prob, initialX);
    }
    else
    if (stdmethodname == "AManPG")
    {
        solver = new AManPG(Prob, initialX);
    }
    else
    {
        std::cout << "Warning: Unrecognized solver: " << stdmethodname << ". Use LRBFGS instead!" << std::endl;
        solver = new LRBFGS(Prob, initialX);
    }
    solver->SetParams(params);

    if(strlen(Handles->isstopped) > 0)
    {
        RJULIA::isstopped = jl_get_function(jl_main_module, Handles->isstopped);
        solver->StopPtr = &RJULIA::juliaInnerStop;
    }
    if(strlen(Handles->LinesearchInput) > 0)
    {
        RJULIA::LinesearchInput = jl_get_function(jl_main_module, Handles->LinesearchInput);
        SolversSMLS *solverSMLS = dynamic_cast<SolversSMLS *> (solver);
        if (solverSMLS != nullptr)
            solverSMLS->LinesearchInput = &RJULIA::juliaLinesearchInput;
    }

    if(Sparams->IsCheckParams != 0)
        solver->CheckParams();
    solver->Run();

    
//    Vector Xopttmp = solver->GetXopt();
//    mexProblem::ObtainMxArrayFromElement(plhs[0], &Xopttmp);
//    plhs[1] = mxCreateDoubleScalar(static_cast<double> (solver->Getfinalfun()));
//
//    SolversSM *solverSM = dynamic_cast<SolversSM *> (solver);
//    if(solverSM != nullptr)
//    {
//        plhs[2] = mxCreateDoubleScalar(static_cast<double> (solverSM->Getnormgf()));
//        plhs[3] = mxCreateDoubleScalar(static_cast<double> (solverSM->Getnormgfgf0()));
//    }
//    SolversNSM *solverNSM = dynamic_cast<SolversNSM *> (solver);
//    if(solverNSM != nullptr)
//    {
//        plhs[2] = mxCreateDoubleScalar(static_cast<double> (solverNSM->Getnormnd()));
//        plhs[3] = mxCreateDoubleScalar(static_cast<double> (solverNSM->Getnormndnd0()));
//    }
//
//    plhs[4] = mxCreateDoubleScalar(static_cast<double> (solver->GetIter()));
//    plhs[5] = mxCreateDoubleScalar(static_cast<double> (solver->Getnf()));
//    plhs[6] = mxCreateDoubleScalar(static_cast<double> (solver->Getng()));
//    plhs[7] = mxCreateDoubleScalar(static_cast<double> (solver->GetnR()));
//    plhs[8] = mxCreateDoubleScalar(static_cast<double> (solver->GetnV()));
//    plhs[9] = mxCreateDoubleScalar(static_cast<double> (solver->GetnVp()));
//    plhs[10] = mxCreateDoubleScalar(static_cast<double> (solver->GetnH()));
//    plhs[11] = mxCreateDoubleScalar(static_cast<double> (solver->GetComTime()));
//    integer lengthSeries = solver->GetlengthSeries();
//    plhs[12] = mxCreateDoubleMatrix(lengthSeries, 1, mxREAL);
//    plhs[13] = mxCreateDoubleMatrix(lengthSeries, 1, mxREAL);
//    plhs[14] = mxCreateDoubleMatrix(lengthSeries, 1, mxREAL);
//
//    double *plhsfun = mxGetPr(plhs[12]), *plhsgrad = mxGetPr(plhs[13]), *plhstime = mxGetPr(plhs[14]); //--, *plhsdist = mxGetPr(plhs[15]);
//    const double *tmpSeries = nullptr;
//    tmpSeries = (solverSM == nullptr) ? solverNSM->GetdirSeries().ObtainReadData() : solverSM->GetgradSeries().ObtainReadData();
//
//    for (integer i = 0; i < lengthSeries; i++)
//    {
//        plhsfun[i] = solver->GetfunSeries().ObtainReadData()[i];
//        plhstime[i] = solver->GettimeSeries().ObtainReadData()[i];
//        plhsgrad[i] = tmpSeries[i];
//    }
//
//    tmp = mexProblem::GetFieldbyName(SolverParams, 0, "IsCheckGradHess");
//    if (tmp != nullptr)
//    {
//        if (fabs(mxGetScalar(tmp)) > std::numeric_limits<realdp>::epsilon()) // if the value is nonzero
//        {
//            Prob->CheckGradHessian(*initialX);
//            Prob->CheckGradHessian(solver->GetXopt());
//        }
//    }
//    plhs[15] = mxCreateDoubleMatrix(4, 1, mxREAL);
//    double *plhseigHess = mxGetPr(plhs[15]);
//    for(integer i = 0; i < 4; i++)
//        plhseigHess[i] = 0;
//
//    if (solver->Verbose >= DETAILED)
//    {
//        Vector MinMaxEigVals1 = Prob->MinMaxEigValHess(*initialX);
//        plhseigHess[0] = MinMaxEigVals1.ObtainReadData()[0];
//        plhseigHess[1] = MinMaxEigVals1.ObtainReadData()[1];
//
//        Vector MinMaxEigVals2 = Prob->MinMaxEigValHess(solver->GetXopt());
//        plhseigHess[2] = MinMaxEigVals2.ObtainReadData()[0];
//        plhseigHess[3] = MinMaxEigVals2.ObtainReadData()[1];
//    }
//
    
    integer lengthSeries = solver->GetlengthSeries();
    integer lengthresult = 1 + initialX->Getlength() + 11 + 3 * lengthSeries + 4;
    double *result = new double[lengthresult];
    result[0] = lengthresult;
    const double *xoptptr = solver->GetXopt().ObtainReadData();
    integer lengthx = initialX->Getlength();
    dcopy_(&lengthx, const_cast<double *>(xoptptr), &GLOBAL::IONE, result + 1, &GLOBAL::IONE);
    result[lengthx + 1] = static_cast<double> (solver->Getfinalfun());

    SolversSM *solverSM = dynamic_cast<SolversSM *> (solver);
    if(solverSM != nullptr)
    {
        result[lengthx + 2] = static_cast<double> (solverSM->Getnormgf());
        result[lengthx + 3] = static_cast<double> (solverSM->Getnormgfgf0());
    }
    SolversNSM *solverNSM = dynamic_cast<SolversNSM *> (solver);
    if(solverNSM != nullptr)
    {
        result[lengthx + 2] = static_cast<double> (solverNSM->Getnormnd());
        result[lengthx + 3] = static_cast<double> (solverNSM->Getnormndnd0());
    }
    result[lengthx + 4] = static_cast<double> (solver->GetIter());
    result[lengthx + 5] = static_cast<double> (solver->Getnf());
    result[lengthx + 6] = static_cast<double> (solver->Getng());
    result[lengthx + 7] = static_cast<double> (solver->GetnR());
    result[lengthx + 8] = static_cast<double> (solver->GetnV());
    result[lengthx + 9] = static_cast<double> (solver->GetnVp());
    result[lengthx + 10] = static_cast<double> (solver->GetnH());
    result[lengthx + 11] = static_cast<double> (solver->GetComTime());
    

    const double *tmpSeries = nullptr;
    tmpSeries = (solverSM == nullptr) ? solverNSM->GetdirSeries().ObtainReadData() : solverSM->GetgradSeries().ObtainReadData();
    
    for (integer i = 0; i < lengthSeries; i++)
        result[lengthx + 12 + i] = solver->GetfunSeries().ObtainReadData()[i];
    
    for (integer i = 0; i < lengthSeries; i++)
        result[lengthx + 12 + lengthSeries + i] = solver->GettimeSeries().ObtainReadData()[i];
    
    for (integer i = 0; i < lengthSeries; i++)
        result[lengthx + 12 + 2 * lengthSeries + i] = tmpSeries[i];
    
    if(Sparams->IsCheckGradHess != 0)
    {
        Prob->CheckGradHessian(*initialX);
        Prob->CheckGradHessian(solver->GetXopt());
    }

    if(Sparams->Verbose >= 3)
    {
        Vector MinMaxEigVals1 = Prob->MinMaxEigValHess(*initialX);
        result[lengthx + 12 + 3 * lengthSeries] = MinMaxEigVals1.ObtainReadData()[0];
        result[lengthx + 12 + 3 * lengthSeries + 1] = MinMaxEigVals1.ObtainReadData()[1];
        
        Vector MinMaxEigVals2 = Prob->MinMaxEigValHess(solver->GetXopt());
        result[lengthx + 12 + 3 * lengthSeries + 2] = MinMaxEigVals2.ObtainReadData()[0];
        result[lengthx + 12 + 3 * lengthSeries + 3] = MinMaxEigVals2.ObtainReadData()[1];
    } else
    {
        result[lengthx + 12 + 3 * lengthSeries] = 0;
        result[lengthx + 12 + 3 * lengthSeries + 1] = 0;
        result[lengthx + 12 + 3 * lengthSeries + 2] = 0;
        result[lengthx + 12 + 3 * lengthSeries + 3] = 0;
    }
    
    delete solver;
    return result;
};

bool ParseManiParams(struct ManiParams *Mparams, Manifold **&manifolds, integer &numoftype, integer *&powsinterval)
//bool ParseManiParams(struct ManiParams *Mparams, Manifold **&manifolds, Element **&elements,
//                     integer *&powsinterval, integer &numoftype, integer &numoftotal)
{
    // Parse ManiParams
    numoftype = Mparams->numoftypes;
    
    powsinterval = new integer[numoftype + 1];
    const char *name = nullptr;
    manifolds = new Manifold *[numoftype];
    integer n, p, m, Params;
    powsinterval[0] = 0;
    
    for (integer i = 0; i < numoftype; i++)
        powsinterval[i + 1] = powsinterval[i] + Mparams->numofmani[i];
    
    PARAMSMAP params;
    for (integer i = 0; i < numoftype; i++)
    {
        name = Mparams->name[i];
        n = Mparams->n[i];
        m = Mparams->m[i];
        p = Mparams->p[i];
        Params = Mparams->paramset[i];
        if(Params != -1)
            params[static_cast<std::string> ("ParamSet")] = Params;

        manifolds[i] = GetAManifold(name, n, m, p);
        manifolds[i]->SetParams(params);

        if (manifolds[i] == nullptr)
        {
            return false;
        }
    }
    return true;
};

Manifold *GetAManifold(const char *name, integer n, integer m, integer p)
{
    if (strcmp(name, "CFixedRankQ2F") == 0)
    {
        return new CFixedRankQ2F(m, n, p);
    }
    else
    if (strcmp(name, "CStiefel") == 0)
    {
        return new CStiefel(n, p);
    }
    else
    if (strcmp(name, "CSymFixedRankQ") == 0)
    {
        return new CSymFixedRankQ(n, p);
    }
    else
    if (strcmp(name, "Euclidean") == 0)
    {
        return new Euclidean(m, n);
    }
    else
    if (strcmp(name, "FixedRankE") == 0)
    {
        return new FixedRankE(m, n, p);
    }
    else
    if (strcmp(name, "FixedRankQ2F") == 0)
    {
        return new FixedRankE(m, n, p);
    }
    else
    if (strcmp(name, "Grassmann") == 0)
    {
        return new Grassmann(n, p);
    }
    else
    if (strcmp(name, "SPDManifold") == 0)
    {
        return new SPDManifold(n);
    }
    else
    if (strcmp(name, "Sphere") == 0)
    {
        return new Sphere(n);
    }
    else
    if (strcmp(name, "Stiefel") == 0)
    {
        return new Stiefel(n, p);
    }
    else
    if (strcmp(name, "SymFixedRankQ") == 0)
    {
        return new SymFixedRankQ(n, p);
    }
    else
    {
        printf("Manifold: %s does not implemented in this library!\n", name);
        return nullptr;
    }
};

#endif // end of DRIVERJULIAPROB
