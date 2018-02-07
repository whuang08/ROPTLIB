function [FinalX, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times, dist] = testLyapunovLR()

    n = 500; k = 2;
    randn('state', 0);%rng(0)

    h = 1/(n+1);
    A = 1/h^2*spdiags(ones(n,1)*[-1 2 -1], [-1 0 1], n, n);
    b = rand(n,1); 

    A = full(A); B = b*b';

    XXex = lyap(A,-B);    
    s = svd(XXex);
    s(1:k+2)   % Xex has very low numerical rank


    % ManiParams.IsCheckParams = 1;
    ManiParams.name = 'LowRank';
    ManiParams.m = n;
    ManiParams.n = n;
    ManiParams.p = k;
    HasHHR = 0;

    problem.M = fixedrankembeddedfactory(n, n, k);
%     X0 = problem.M.rand();  % very very bad performance
    % For checking, let us try the best rank k approx of Xex and perturb a
    % little
    [UU,SS,VV] =  svd(XXex);
    X0.U = UU(:,1:k); X0.V = VV(:,1:k); X0.S = SS(1:k,1:k);
    Noise = problem.M.randvec(X0);  sz = 1e-1;
    X0 = problem.M.retr(X0,Noise,sz);
    initialX.main = [X0.U(:); X0.S(:); X0.V(:)];

    fhandle = @(x)f(x, A, B, n, k);
    gfhandle = @(x)gf(x, A, B, n, k);
    Hesshandle = @(x, eta)Hess(x, eta, A, B, n, k, problem);

%     SolverParams.method = 'LRBFGS';
    SolverParams.method = 'LRTRSR1';
%     SolverParams.method = 'RTRNewton';
%     SolverParams.method = 'RCG';
    % SolverParams.IsCheckParams = 1;
    SolverParams.DEBUG = 2;
    % SolverParams.Min_Iteration = 100;
%     SolverParams.IsCheckGradHess = 1;
    SolverParams.Max_Iteration = 1000;
    SolverParams.OutputGap = 1;
    
    [FinalX, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times] = DriverOPT(fhandle, gfhandle, Hesshandle, SolverParams, ManiParams, HasHHR, initialX);
%     [FinalX, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times, dists] = DriverOPT(fhandle, gfhandle, Hesshandle, SolverParams, ManiParams, HasHHR, initialX, FinalX);

    p = k;
    U = reshape(FinalX.main(1 : n * p), n, p);
    S = reshape(FinalX.main(n * p + 1 : n * p + p * p), p, p);
    V = reshape(FinalX.main(n * p + p * p + 1 : n * p + p * p + n * p), n, p);
    XX = U * S * V';
    error = norm( XX - XXex ) / norm( XXex )    

end

function [output, x] = f(x, A, B, n, p)
    %if ~isfield(store, 'X_on_omega')
    %    store.X_on_omega = partXY(X.S'*X.U',X.V', prob.Omega_i, prob.Omega_j, prob.m)'; 
    %end
    %X_on_omega = store.X_on_omega;

    U = reshape(x.main(1 : n * p), n, p);
    S = reshape(x.main(n * p + 1 : n * p + p * p), p, p);
    V = reshape(x.main(n * p + p * p + 1 : n * p + p * p + n * p), n, p);
    XX = U * S * V';
    output = 0.5*trace(XX'*(A*XX+XX*A-2*B));  
%     outputinf = output
    
end

function [output, x] = gf(x, A, B, n, p)
    %         if ~isfield(store, 'X_on_omega')
    %             store.X_on_omega = partXY(X.S'*X.U',X.V', prob.Omega_i, prob.Omega_j, prob.m)'; 
    %         end
    %         X_on_omega = store.X_on_omega;
    U = reshape(x.main(1 : n * p), n, p);
    S = reshape(x.main(n * p + 1 : n * p + p * p), p, p);
    V = reshape(x.main(n * p + p * p + 1 : n * p + p * p + n * p), n, p);
    XX = U * S * V';
    G = A*XX+XX*A-B;
    output.EucRep = [0; G(:)];
%     'gf'
%     dS = U'* G * V;
%     dU = (G * V - U * (U' * (G * V))) / S;
%     dV = S \ (U' * G - ((U' * G) * V) * V');
%     output.main = [dU(:); dS(:); dV(:)];
end

function [output, x] = Hess(x, eta, A, B, n, p, problem)
    X.U = reshape(x.main(1 : n * p), n, p);
    X.S = reshape(x.main(n * p + 1 : n * p + p * p), p, p);
    X.V = reshape(x.main(n * p + p * p + 1 : n * p + p * p + n * p), n, p);
    dU = reshape(eta.main(1 : n * p), n, p);
    dS = reshape(eta.main(n * p + 1 : n * p + p * p), p, p);
    dV = reshape(eta.main(n * p + p * p + 1 : n * p + p * p + n * p), n, p);
    H.Up = dU * X.S;
    H.Vp = dV * X.S';
    H.M = dS;
    
    ambient_H = problem.M.tangent2ambient(X, H);
    Xdot = ambient_H.U*ambient_H.S*ambient_H.V';        
    ehess = A*Xdot+Xdot*A;
    output.EucRep = [0; ehess(:)];
end
