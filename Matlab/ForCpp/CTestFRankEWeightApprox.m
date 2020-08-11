function MTestWeightedLowRank()

    m = 4; n = 4; r = 2;
    A = randn(m, n);
    O = orth(randn(m * n, m * n));
    D = diag(rand(m * n, 1) + 0.1);
    W = O * D * O';
    U = orth(randn(m, r));
    D = randn(r, r);
    V = orth(randn(n, r));
    Xinitial.main = U * D * V';
    Xinitial.U = U;
    Xinitial.D = D;
    Xinitial.V = V;
    HasHHR = 0;
    
    SolverParams.method = 'LRBFGS';
%     SolverParams.method = 'LRTRSR1';
%     SolverParams.method = 'RTRNewton';
%     SolverParams.method = 'RCG';
    SolverParams.IsCheckParams = 1;
    SolverParams.Max_Iteration = 500;
    SolverParams.LengthSY = 4;
    SolverParams.Verbose = 2;
    SolverParams.LMrestart = 0;
    SolverParams.OutputGap = 10;
%     SolverParams.Accuracy = 1e-6;
%     SolverParams.
%     SolverParams.Tolerance = 1e-9;% * norm(B);
    SolverParams.Finalstepsize = 1;
    SolverParams.Num_pre_funs = 2;
    SolverParams.PreFunsAccuracy = 1e6;
%     SolverParams.IsCheckGradHess = 1;

    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times] = ...
        TestFRankEWeightApprox(A, W, Xinitial, HasHHR, SolverParams);
end
