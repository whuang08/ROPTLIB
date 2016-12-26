function MTestEucPosSpCd()
    r = floor(rand() * 10000000)
%     r = 1;
    r = 9898359
%     r = 7740966
    rand('state', r);
    randn('state', r);
    dim = 10;
    num = 200;
    N = 1;
    
    alphatrue = max(rand(num, N) - 0.2, 0);


    Bs = zeros(dim, dim, num);
    for i = 1 : num
        tmp = randn(dim, dim * 10);
        Bs(:, :, i) = tmp * tmp';
    end
    Balpha = reshape(reshape(Bs, dim * dim, num) * alphatrue, dim, dim, N);
    
    Ls = zeros(dim, dim, N);
    for i = 1 : N
        Ls(:, :, i) = chol(Balpha(:, :, i), 'lower');
    end
    
    pinvBs = pinv(reshape(Bs, dim * dim, []));
    alphainit = rand(num, N);
    for i = 1 : N
%         tmpv = pinvBs * reshape(Balpha(:, :, i), [], 1);
%         alphainit(:, i) = max(tmpv, 0);
%         alphainit(:, i) = max(tmpv, 0) / sum(abs(tmpv)) * sum(max(tmpv, 0));

        a = norm(pinvBs * reshape(Balpha(:, :, i), [], 1))
        alphainit(:, i) = ones(num, 1) / sqrt(num) * a;
%         alphainit(:, i) = alphainit(:, i) / norm(alphainit(:, i)) * a;
    end
    
    SolverParams.method = 'LRBFGS';
    SolverParams.IsCheckParams = 1;
    SolverParams.Max_Iteration = 500;
    SolverParams.OutputGap = 100;
    SolverParams.LengthSY = 4;
    SolverParams.Num_pre_funs = 0;
    SolverParams.InitSteptype = 3;
    SolverParams.DEBUG = 2;
    SolverParams.Stop_Criterion = 2;
%     SolverParams.IsCheckGradHess = 1;
    HasHHR = 0;
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times] = TestEucPosSpCd(Ls, Bs, alphainit, 0.1, HasHHR, SolverParams);
    
%     norm(reshape(Xopt.main, num, N) - alphatrue) / norm(alphatrue)
%     reshape(Xopt.main, num, N)
    
    
    figure(1);clf
    semilogy(1:length(funs), funs, 'ob-');
    hold on
    figure(2);clf
    semilogy(1:length(grads), grads, 'ob-');
    hold on
    
    SolverParams.method = 'RSD';
    SolverParams.IsCheckParams = 1;
    SolverParams.Max_Iteration = 500;
    SolverParams.OutputGap = 100;
    SolverParams.LengthSY = 4;
    SolverParams.Num_pre_funs = 0;
    SolverParams.InitSteptype = 1;
    SolverParams.DEBUG = 2;
    SolverParams.Stop_Criterion = 2;
%     SolverParams.IsCheckGradHess = 1;
    HasHHR = 0;
    [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times] = TestEucPosSpCd(Ls, Bs, alphainit, 0.1, HasHHR, SolverParams);

    figure(1);
    semilogy(1:length(funs), funs, 'xr-');
    hold on
    legend('LRBFGS', 'RSD');
    title('f');
    figure(2);
    semilogy(1:length(grads), grads, 'xr-');
    hold on
    legend('LRBFGS', 'RSD');
    title('gf');
end
