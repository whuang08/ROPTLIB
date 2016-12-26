function MTestStieSparseBrockett()
    r = floor(rand() * 100000);
    r = 2
    fprintf('seed:%d\n', r);
    rand('state', r);
    randn('state', r);
    n = 6;
    p = 2;
%     B = sprandn(n, n, 1/n);
%     B = B + B';
    B = sparse(n, n);
    for i = 1 : n
        B(i, i) = i;
    end
%     [V, D] = eigs(B, p, 'SA');
%     tmp = randn(n, p); tmp = tmp / norm(tmp);
%     Xinitial = orth(V + tmp * 0.001);
    Xinitial = orth(randn(n, p));
    D = ones(p, 1);%(p:-1:1)';%
    
    SolverParams.method = 'RSD';
%     SolverParams.method = 'RTRSR1';
%     SolverParams.method = 'RTRNewton';
    SolverParams.IsCheckParams = 1;
    SolverParams.Max_Iteration = 2;
    SolverParams.LengthSY = 4;
    SolverParams.DEBUG = 2;
    SolverParams.OutputGap = 1;
    SolverParams.LS_ratio1 = 0.5;
    SolverParams.LS_ratio2 = 0.5;
    SolverParams.InitSteptype = 0;
    SolverParams.InitStepsize = 1;
    
    SolverParams.LineSearch_LS = 4;
    SolverParams.LinesearchInput = @(x, eta, t0, s0)LinesearchInput(x, eta, t0, s0, B, D, n, p);
    SolverParams.Accuracy = 1e-4;
    SolverParams.Tolerance = 1e-4;% * norm(B);
    SolverParams.Finalstepsize = 1;
    HasHHR = 0;
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime] = TestStieSparseBrockett(B, D, Xinitial, HasHHR, 1, SolverParams);
%     [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime] = TestStieBrockett(full(B), D, Xinitial, HasHHR, 1, SolverParams);
    HasHHR = 0;
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime] = TestStieSparseBrockett(B, D, Xinitial, HasHHR, 3, SolverParams);
%     [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime] = TestStieBrockett(full(B), D, Xinitial, HasHHR, 3, SolverParams);
end

function output = LinesearchInput(x, eta, t0, s0, B, D, n, p)
    tmp = B * reshape(eta.main, n, p) * diag(D);
    fprintf('Tr(x B eta D)%f, Tr(eta B eta D) %f\n', sum(tmp(:)'*x.main), sum(tmp(:)'*eta.main));%%---
    output = - sum(tmp(:)'*x.main) / sum(tmp(:)'*eta.main); %fixing the stepsize to be one does not guarantee convergence.
    xx = reshape(x.main, n, p)
    eet = reshape(output * eta.main, n, p)
    
    2 * B * xx - 2 * xx * (xx' * B * xx)
    
end
