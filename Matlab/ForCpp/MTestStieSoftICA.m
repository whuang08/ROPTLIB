function MTestStieSoftICA()

    rand('state',1);
    randn('state',1);
%     n = 1000; p = 3; N = 5;
%     Cs = randn(n, n, N);
%     for i = 1 : N
%         Cs(:, :, i) = Cs(:, :, i) + Cs(:, :, i)';
%     end
%     X = orth(randn(n, p));
%     SolverParams.method = 'LRBFGS';
%     SolverParams.DEBUG = 1;
%     HasHHR = 0;
%     [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funSeries, gradSeries, timeSeries] = TestStieSoftICA(Cs, X, n, p, HasHHR, 1, SolverParams);

    n = 10; p = 3; N = 5;
    Cs = rand(n, n, N);
    for i = 1 : N
        Cs(:, :, i) = Cs(:, :, i) + Cs(:, :, i)' + diag(1:n);
    end
    X = orth(rand(n, p));
    SolverParams.method = 'LRBFGS';
    SolverParams.DEBUG = 1;
    SolverParams.Max_Iteration = 2000;
    HasHHR = 0;
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funSeries, gradSeries, timeSeries] = TestStieSoftICA(Cs, X, n, p, HasHHR, 2, SolverParams);
% 
%     SolverParams.method = 'RCG';
%     SolverParams.DEBUG = 1;
%     SolverParams.Num_pre_funs = 2;
%     SolverParams.Max_Iteration = 2000;
%     SolverParams.RCGmethod = 4;
%     HasHHR = 0;
%     [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funSeries, gradSeries, timeSeries] = TestStieSoftICA(Cs, X, n, p, HasHHR, 1, SolverParams);
end
