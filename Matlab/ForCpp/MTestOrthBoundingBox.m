function MTestOrthBoundingBox()
    r = floor(rand() * 100000);
%     r = 1
    r = 81472
    fprintf('seed:%d\n', r);
    rand('state', r);
    randn('state', r);
    d = 3;
    n = 100;
    Xinitial = orth(randn(d, d));
    E = rand(d, n);
    
    SolverParams.method = 'RBFGSLPSub';
    SolverParams.lambdaLower = 1e-2;
    SolverParams.lambdaUpper = 1e2;
    SolverParams.IsCheckParams = 1;
    SolverParams.Max_Iteration = 50;
    SolverParams.NumExtraGF = floor(d * (d - 1) / 2);
    SolverParams.LengthSY = 0;
    SolverParams.DEBUG = 2;
    SolverParams.ParamSet = 2;
%     SolverParams.Eps = 1e-4;
%     SolverParams.Del = 1e-4;
%     SolverParams.Tolerance = 1e-10;
    SolverParams.Stop_Criterion = 3;
    SolverParams.LineSearch_LS = 4;
    SolverParams.OutputGap = 10;
%     SolverParams.Diffx = 1e-7;
    SolverParams.Finalstepsize = 1;
    HasHHR = 1;
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funSeries, gradSeries, timeSeries] = TestOrthBoundingBox(E, Xinitial, HasHHR, 1, SolverParams);
%     [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funSeries, gradSeries, timeSeries] = TestOrthBoundingBoxRMethod(E, Xinitial, HasHHR, 1, SolverParams);
%     funSeries
%     xopt = reshape(Xopt.main, d, d);
%     OE = xopt * E;
%     scatter3(OE(1, :), OE(2, :), OE(3, :), '.');
end
