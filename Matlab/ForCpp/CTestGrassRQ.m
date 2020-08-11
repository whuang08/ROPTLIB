function MTestGrassRQ()
    seed = floor(rand() * 100000);
    seed = 2;
    fprintf('MTestGrassRQ seed:%d\n', seed);
    rand('state', seed);
    randn('state', seed);
    n = 100;
    p = 4;
    B = randn(n, n);
    B = B + B';
    Xinitial = orth(randn(n, p));
    SolverParams.method = 'LRBFGS';
%     SolverParams.method = 'RTRSR1';
%     SolverParams.method = 'RTRNewton';
    SolverParams.IsCheckParams = 1;
    SolverParams.Max_Iteration = 500;
    SolverParams.LengthSY = 4;
    SolverParams.Verbose = 1;
    SolverParams.Accuracy = 1e-6;
%     SolverParams.Tolerance = 1e-8;%1e-10;% * norm(B);
%     SolverParams.Finalstepsize = 1;
%     SolverParams.InitSteptype = 3;
%     SolverParams.IsCheckGradHess = 1;
    HasHHR = 0;
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime] = TestGrassRQ(B, Xinitial, HasHHR, SolverParams);
end
