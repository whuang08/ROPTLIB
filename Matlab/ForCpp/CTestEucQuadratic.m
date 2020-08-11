function MTestEucQuadratic()
    seed = floor(rand() * 100000);
%     seed = 3;
    fprintf('MTestEucQuadratic seed:%d\n', seed);
    rand('state', seed);
    randn('state', seed);
    dim = 100;
    O = orth(randn(100, 100));
    D = diag([0.001 * (rand(50, 1) + 0.5); 1000 * (rand(50, 1) + 0.5)]);
    M = O' * D * O;
    Xinitial = randn(dim, 1);
    
%     SolverParams.method = 'RBFGS';
%     SolverParams.method = 'RBroydenFamily';
    SolverParams.method = 'RTRNewton';
%     SolverParams.method = 'RCG';
%     SolverParams.method = 'LRBFGS';
    SolverParams.IsCheckParams = 1;
    SolverParams.Max_Iteration = 2000;
    SolverParams.LengthSY = 4;
    SolverParams.Verbose = 2;
    SolverParams.LMrestart = 0;
    SolverParams.Accuracy = 1e-6;
    SolverParams.Tolerance = 1e-9;
    SolverParams.OutputGap = 100;
    SolverParams.Finalstepsize = 1;
    SolverParams.Num_pre_funs = 1;
    SolverParams.PreFunsAccuracy = 1e6;
    HasHHR = 0;
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times] = TestEucQuadratic(M, Xinitial, HasHHR, SolverParams);

end
