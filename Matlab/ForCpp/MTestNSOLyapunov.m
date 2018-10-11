function MTestNSOLyapunov()
    seed = floor(rand() * 100000);
    seed = 2;
    fprintf('MTestNSOLyapunov seed:%d\n', seed);
    rand('state', seed);
    randn('state', seed);
%     n = 50;
%     p = 3;
%     A = randn(n, n); A = A' * A;
%     O = orth(randn(n, n));D=diag(logspace(1,2,n));M=O*D*O';
%     c = randn(n, 1); C = c * c';
    
    n = 500; p = 2;
    randn('state', 0)
    rand('state', 0)
    h = 1/(n+1);
    A = 1/h^2*spdiags(ones(n,1)*[-1 2 -1], [-1 0 1], n, n);
    A = full(A); 
    M = eye(n);
    b = rand(n,1); C = b*b';
    XXex = lyap(A,-C);    
    s = svd(XXex);
    s(1:p+2)   % Xex has very low numerical rank
    
    
    
    Xinitial = randn(n, p);
    metric = 1;
    
%     SolverParams.method = 'RSD';
    SolverParams.method = 'LRBFGS';
%      SolverParams.method = 'RCG';
%     SolverParams.method = 'LRTRSR1';
%     SolverParams.method = 'RTRSR1';
%     SolverParams.method = 'RTRNewton';
%     SolverParams.IsCheckGradHess = 1;
    SolverParams.IsCheckParams = 1;
    SolverParams.Max_Iteration = 4000;
    SolverParams.OutputGap = 50;
    SolverParams.LengthSY = 4;
    SolverParams.DEBUG = 2;
%     SolverParams.InitSteptype = 1;
    SolverParams.Tolerance = 1e-7;
    HasHHR = 0;
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times, dists] = TestNSOLyapunov(A, M, C, Xinitial, HasHHR, metric, SolverParams);
end
