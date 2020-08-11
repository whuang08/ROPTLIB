function MTestSphereSparestVector()
    seed = floor(rand() * 100000);
    seed = 81472;
    fprintf('MTestSphereSparestVector seed:%d\n', seed);
    rand('state', seed);
    randn('state', seed);
    m = 100;
    n = 10;
    Xinitial = orth(randn(n, 1));
    Q = randn(m, n);
    Q(:, 1) = zeros(m, 1);
    Q(1, 1) = 1;
    
    SolverParams.method = 'RGS';
    SolverParams.lambdaLower = 1e-2;
    SolverParams.lambdaUpper = 1e2;
    SolverParams.IsCheckParams = 0;
    SolverParams.Max_Iteration = 200;
%     SolverParams.NumExtraGF = 2;
    SolverParams.LengthSY = 0;
    SolverParams.Verbose = 1;
    SolverParams.ParamSet = 1;
%     SolverParams.Eps = 1e-4;
%     SolverParams.Del = 1e-4;
%     SolverParams.Tolerance = 1e-10;
    SolverParams.Stop_Criterion = 0;
%     SolverParams.LineSearch_LS = 4;
    SolverParams.OutputGap = 1;
%     SolverParams.Diffx = 1e-7;
%     SolverParams.Finalstepsize = 1;
    SolverParams.IsCheckParams = 1;
    SolverParams.Tolerance = 1e-6;
    HasHHR = 1;
    ParamSet = 1;
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funSeries, gradSeries, timeSeries] = TestSphereSparsestVector(Q, Xinitial, HasHHR, ParamSet, SolverParams);
%     sum(abs(Q * Xopt.main)<1e-5)

%     funSeries
%     xopt = reshape(Xopt.main, d, d);
%     OE = xopt * E;
%     scatter3(OE(1, :), OE(2, :), OE(3, :), '.');
end
