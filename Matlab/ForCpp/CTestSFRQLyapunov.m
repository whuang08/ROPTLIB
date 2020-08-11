function MTestNSOLyapunov()

    n = 100; k = 2;
    seed = floor(rand() * 10000000)
%     seed = 5877728
    randn('state', seed);%rng(0)
    rand('state', seed);%rng(0)

%     A = randn(n, n);A=A*A';
    h = 1/(n+1);
    A = spdiags(ones(n,1)*[-1 2 -1], [-1 0 1], n, n); %(1/h^2)*
%     M = spdiags((1 + ((1:n)' / n) * 300), [0], n, n);
    M = spdiags(ones(n, 1), [0], n, n);
%     A = A / norm(A, 'fro') * n;
    b = randn(n,1); 
    b = b;% / norm(b);
    C = b;

%     XXex = lyap(A,-C * C');    
%     s = svd(XXex);
%     s(1:p+2)   % Xex has very low numerical rank
%     
    
    Xinitial = randn(n, k);
    ParamSet = 2;
    
%     SolverParams.method = 'RSD';
    SolverParams.method = 'LRBFGS';
%      SolverParams.method = 'RCG';
%     SolverParams.method = 'LRTRSR1';
%     SolverParams.method = 'RTRSR1';
%     SolverParams.method = 'RTRNewton';
%     SolverParams.IsCheckGradHess = 1;
    SolverParams.IsCheckParams = 0;
    SolverParams.Max_Iteration = 5000;
    SolverParams.OutputGap = 10;
%     SolverParams.LengthSY = 0;
%     SolverParams.BBratio = -1;%initial Hessian = Identity
    SolverParams.Verbose = 3;
%     SolverParams.InitSteptype = 0;
%     SolverParams.LineSearch_LS = 0;
%     SolverParams.Initstepsize = 0.001;
%     SolverParams.Tolerance = 1e-8;
    HasHHR = 0;
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times, Heigs] = TestSFRQLyapunov(A, M, C, Xinitial, HasHHR, ParamSet, SolverParams);
    Heigs
    %     svd(reshape(Xopt.main, n, k)).^2
%     solverTmp.acceptedstepsize
end
