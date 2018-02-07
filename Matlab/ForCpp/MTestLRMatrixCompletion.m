function MTestLRMatrixCompletion()
    seed = floor(rand() * 100000);
    seed = 2;
    fprintf('MTestLRMatrixCompletion seed:%d\n', seed);
    rand('state', seed);
    randn('state', seed);
    m = 100;
    n = 100;
    r = 5;
    
    G = randn(m, r);
    H = randn(n, r);
    B = G * H';
    OS = 3;
%     [(m + n - r) * r * OS, m * n]
    nz = min((m + n - r) * r * OS, m * n);
    
    vidx = randperm(m * n, nz);
    [ir, jc] = ind2sub([m, n], vidx);
    A = sparse(ir, jc, B(vidx));

    [U, D, V] = svds(full(A),r);
    Xinitial = [U(:); D(:); V(:)];
    
    SolverParams.method = 'LRBFGS';
%     SolverParams.method = 'RTRSR1';
%     SolverParams.method = 'RTRNewton';
    SolverParams.IsCheckParams = 0;
    SolverParams.Max_Iteration = 1000;
    SolverParams.OutputGap = 100;
    SolverParams.LengthSY = 4;
    SolverParams.DEBUG = 1;
    HasHHR = 0;
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime] = TestLRMatrixCompletion(A, Xinitial, r, HasHHR, SolverParams);
end
