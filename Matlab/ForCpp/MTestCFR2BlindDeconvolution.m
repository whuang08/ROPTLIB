function MTestLRMatrixCompletion()
    seed = floor(rand() * 100000);
    seed = 2;
    fprintf('seed:%d\n', seed);
    rand('state', seed);
    randn('state', seed);
    K = 100;
    N = 100;
    L = round(2*(K+N)); 
    r = 1;
    
    idx2 = randperm(L);
    idx2 = idx2(1:K);
    Breal = speye(L);
    B = Breal(:,idx2);
    
%     B = randn(2 * L, K);
    C = complex(randn(L, N), randn(L, N));
    y = complex(randn(L, 1), randn(L, 1));
    
    Xinitial = randn(2*(K + N) * r);
    
    SolverParams.method = 'LRBFGS';
%     SolverParams.method = 'RTRSR1';
%     SolverParams.method = 'RTRNewton';
    SolverParams.IsCheckParams = 1;
    SolverParams.Max_Iteration = 1000;
    SolverParams.OutputGap = 100;
    SolverParams.LengthSY = 4;
    SolverParams.DEBUG = 2;
    HasHHR = 0;
    rho = 0.;
    d = 1;
    mu = 1;
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times, dists] = TestCFR2BlindDeconvolution(y, B, C, Xinitial, r, HasHHR, SolverParams, rho, d, mu);
end
