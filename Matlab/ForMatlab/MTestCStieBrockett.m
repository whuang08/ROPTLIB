clear;
n = 5;
p = 2;
randn('state', 1);
rand('state', 1);
B = complex(randn(n, n), randn(n, n)); B = B + B';
D = sparse(diag(p : -1 : 1));
fhandle = @(x)f(x, B, D);
gfhandle = @(x)gf(x, B, D);
Hesshandle = @(x, eta)Hess(x, eta, B, D);
PreConhandle = @(x, eta)PreCon(x, eta);

% SolverParams.method = 'RTRNewton';
SolverParams.method = 'LRBFGS';
SolverParams.IsCheckParams = 1;
SolverParams.Verbose = 1;
% SolverParams.Min_Iteration = 100;
SolverParams.IsCheckGradHess = 0;
SolverParams.Tolerance = 1e-10;
SolverParams.Accuracy = 1e-6;
SolverParams.Max_Iteration = 100;
SolverParams.OutputGap = 1;
% SolverParams.InitSteptype = 1;
SolverParams.IsStopped = @IsStopped;

SolverParams.LineSearch_LS = 0;
SolverParams.LinesearchInput = @LinesearchInput;

ManiParams.IsCheckParams = 1;
ManiParams.name = 'CStiefel';
ManiParams.n = n;
ManiParams.p = p;
ManiParams.ParamSet = 4;
ManiParams.IsCheckParams = 1;
HasHHR = 0;

initialX.main = orth(complex(randn(n, p), randn(n, p)));
% use the function handles
[FinalX, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times, eigHess] = DriverOPT(fhandle, gfhandle, Hesshandle, PreConhandle, SolverParams, ManiParams, HasHHR, initialX);
% use numerical gradient and numerical Hessians
[FinalX, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times, eigHess] = DriverOPT(fhandle, [], [], [], SolverParams, ManiParams, HasHHR, initialX);

function output = LinesearchInput(x, eta, t0, s0)
    output = 1;
end

function output = IsStopped(x, funSeries, ngf, ngf0)
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

function [output, x] = PreCon(x, eta)
    output = eta;
end
