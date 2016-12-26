function [FinalX, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times, dist] = testSimpleExample()
    n = 10;
    p = 4;
    B = randn(n, n); B = B + B';
    D = sparse(diag(p : -1 : 1));
    fhandle = @(x)f(x, B, D);
    gfhandle = @(x)gf(x, B, D);
    Hesshandle = @(x, eta)Hess(x, eta, B, D);

    SolverParams.method = 'LRBFGS';
    % SolverParams.IsCheckParams = 1;
    SolverParams.DEBUG = 2;
    % SolverParams.Min_Iteration = 100;
    % SolverParams.IsCheckGradHess = 1;
    SolverParams.Max_Iteration = 500;
    SolverParams.OutputGap = 30;
    SolverParams.IsStopped = @IsStopped;
    
    SolverParams.LineSearch_LS = 0;
    SolverParams.LinesearchInput = @LinesearchInput;

    % ManiParams.IsCheckParams = 1;
    ManiParams.name = 'Stiefel';
    ManiParams.n = n;
    ManiParams.p = p;
    ManiParams.ParamSet = 1;
    HasHHR = 0;

    initialX.main = orth(randn(n, p));
    [FinalX, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times] = DriverOPT(fhandle, gfhandle, Hesshandle, SolverParams, ManiParams, HasHHR, initialX);
    [FinalX, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times, dists] = DriverOPT(fhandle, gfhandle, Hesshandle, SolverParams, ManiParams, HasHHR, initialX, FinalX);
end

function output = LinesearchInput(x, eta, t0, s0)
    output = 1;
end

function output = IsStopped(x, gf, f, ngf, ngf0)
    output = ngf / ngf0 < 1e-5;
end

function [output, x] = f(x, B, D)
x.BUD = B * x.main * D;
output = x.main(:)' * x.BUD(:);
end

function [output, x] = gf(x, B, D)
output.main = 2 * x.BUD;
end

function [output, x] = Hess(x, eta, B, D)
output.main = 2 * B * eta.main * D;
end
