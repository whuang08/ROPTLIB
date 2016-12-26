function MTestGrassRQ()
    r = rand() * 100000;
%     r = 1
    fprintf('seed:%d\n', r);
    rand('state', r);
    randn('state', r);
    n = 1000;
    p = 4;
    B = randn(n, n);
    B = B + B';
    Xinitial = orth(randn(n, p));
    SolverParams.method = 'LRBFGS';
%     SolverParams.method = 'RTRSR1';
%     SolverParams.method = 'RTRNewton';
    SolverParams.IsCheckParams = 1;
    SolverParams.Max_Iteration = 500;
    SolverParams.LengthSY = 4;
    SolverParams.DEBUG = 1;
    SolverParams.Accuracy = 1e-6;
    SolverParams.Tolerance = 1e-8;%1e-10;% * norm(B);
    SolverParams.Finalstepsize = 1;
    SolverParams.InitSteptype = 3;
    HasHHR = 0;
    tic
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime] = TestGrassRQ(B, Xinitial, HasHHR, SolverParams);
    toc
    Xopt.main = reshape(Xopt.main, size(Xinitial));
%     (Xopt.main' * B * Xopt.main)
%     (Xopt.main' * Xopt.main)
%     norm(B * Xopt.main - Xopt.main * (Xopt.main' * B * Xopt.main), 'fro')
%     for i = 1 : p
%         tmp = B * Xopt.main(:, i);
%         lambdai = - norm(tmp);
%         error(i) = norm(tmp - lambdai * Xopt.main(:, i));
%     end
%     error
%     norm(error) / norm(B)
%     opts.issym = 1;
%     opts.isreal = 1;
%     opts.tol = 1e-10;
%     opts.disp = 0;
%     tic
%     [V, D] = eigs(B, p, 'SA', opts);
%     toc
%     (sum(diag(D)) - f)/f
%     (sum(diag(D)) - f)
end
