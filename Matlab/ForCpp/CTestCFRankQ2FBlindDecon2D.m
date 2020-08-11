% function MTestCFR2BlindDecon2D()
    seed = floor(rand() * 100000);
    seed = 2;
    fprintf('MTestCFR2BlindDecon2D seed:%d\n', seed);
    rand('state', seed);
    randn('state', seed);
    n1 = 32;
    n2 = 32;
    L = n1 * n2;
    r = 1;
    
%     idx2 = randperm(L);
%     idx2 = idx2(1:L);
    B = spdiags(complex(ones(L, 1), eps * ones(L, 1)), 0, L, L);
%     B = Breal(:,idx2);
    
    C = spdiags(complex(ones(L, 1), eps * ones(L, 1)), 0, L, L);
    y = complex(randn(L, 1), randn(L, 1));
    
    Xinitial = randn(2*(L + L) * r);
    
%     SolverParams.method = 'RSD';
%     SolverParams.method = 'LRBFGS';
%     SolverParams.method = 'RTRSR1';
    SolverParams.method = 'RTRNewton';
    SolverParams.IsCheckParams = 0;
%     SolverParams.IsCheckGradHess = 1;
    SolverParams.Max_Iteration = 100;
    SolverParams.OutputGap = 1;
    SolverParams.LengthSY = 4;
    SolverParams.Verbose = 2;
%     SolverParams.InitSteptype = 1;
    HasHHR = 0;
    rho = 0.;
    d = 1;
    mu = 1;
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times, dists] = TestCFRankQ2FBlindDecon2D(y, B, C, Xinitial, n1, n2, r, HasHHR, SolverParams, rho, d, mu);
% end
