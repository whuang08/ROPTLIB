function MTestStieBrockett()
    seed = floor(rand() * 100000);
    seed = 1;
    fprintf('MTestStieBrockett seed:%d\n', seed);
    rand('state', seed);
    randn('state', seed);
    n = 10;
    p = 4;
    m = 4;
    q = 2;
    
    B1 = randn(n, n);
    B1 = B1 + B1';
    B2 = randn(n, n);
    B2 = B2 + B2';
    B3 = randn(m, m);
    B3 = B3 + B3';
    
    initX1 = orth(randn(n, p));
    initX2 = orth(randn(n, p));
    initX3 = orth(randn(m, q));
    initX = [initX1(:); initX2(:); initX3(:)];
    
    D1 = (p:-1:1)';
    D2 = (p:-1:1)';
    D3 = (q:-1:1)';
    
    SolverParams.method = 'LRBFGS';
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
    HasHHR = 0;
    paramset = 1;
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times] = ...
        TestProdStieSumBrockett(B1, D1, B2, D2, B3, D3, initX, HasHHR, paramset, SolverParams);
end
