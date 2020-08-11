% function MTestStieSPCA()
    seed = floor(rand() * 100000);
    seed = 20;
    fprintf('MTestFRankESparseApprox seed:%d\n', seed);
    rng(seed);
    
    n = 20;
    m = 20;
    r = 5;
    A = randn(m, n);
    U = orth(randn(m, r));
    D = randn(r, r);
    V = orth(randn(n, r));
    Xinitial.main = U * D * V';
    Xinitial.U = U;
    Xinitial.D = D;
    Xinitial.V = V;
    lambda = 1;
    
    SolverParams.method = 'AManPG';
    SolverParams.IsCheckParams = 1;
    SolverParams.RPGVariant = 0; %0: RPG without adaptive stepsize, 1: RPG with adaptive stepsize
    SolverParams.LengthW = 1;
    SolverParams.OutputGap = 20;
    SolverParams.Max_Iteration = 1000;
    SolverParams.DEBUG = 1;
    [xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times] = ...
        TestFRankESparseApprox(A, lambda, Xinitial, SolverParams);
    
% end
