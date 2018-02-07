function MTestStieBrockett()
    seed = floor(rand() * 100000);
    seed = 1;
    fprintf('MTestStieBrockett seed:%d\n', seed);
    rand('state', seed);
    randn('state', seed);
    n = 100;
    p = 4;
    Xinitial = orth(randn(n, p));
    B = randn(n, n);
    B = B + B';
    
    D = ones(p, 1);%(p:-1:1)'; %
    SolverParams.method = 'LRBFGS';
%     SolverParams.method = 'RTRNewton';
%     SolverParams.method = 'RCG';
    SolverParams.IsCheckParams = 0;
    SolverParams.Max_Iteration = 2000;
    SolverParams.LengthSY = 4;
    SolverParams.DEBUG = 1;
    SolverParams.Accuracy = 1e-6;
    SolverParams.Tolerance = 1e-9;% * norm(B);
    SolverParams.Finalstepsize = 1;
    HasHHR = 0;
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times] = TestStieBrockett(B, D, Xinitial, HasHHR, 1, SolverParams);
%     SolverParams.OutputGap = 40;
%     SolverParams.DEBUG = 2;
%     [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times, dists] = TestStieBrockett(B, D, Xinitial, HasHHR, 1, SolverParams, Xopt.main);
end
