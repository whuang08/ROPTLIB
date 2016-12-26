function MTestSPDTensorDL()
    r = floor(rand() * 10000000)
%     r = 1
%     r = 7740966
    rand('state', r);
    randn('state', r);
    dim = 3;%% 20;
    num = 100;
    N = 500;
    numatoms = 10;%floor(num * 1);
    
    % True dictionary: Xtrue
    Xtrue = zeros(dim, dim, num);
    for i = 1 : num
        tmp = randn(dim, dim * 10);
        Xtrue(:, :, i) = tmp * tmp';
    end
    % True sparse codes: alpha
    alpha = zeros(num, N);
    for i = 1 : N
        values = rand(numatoms, 1);
        idx = randperm(num, numatoms);
        for j = 1 : numatoms
            alpha(idx(j), i) = values(j);
        end
    end
    
    % The input data tensor
    Xalpha = reshape(Xtrue, dim * dim, num) * alpha;
    As = reshape(Xalpha, dim, dim, N);
    Ls = zeros(dim, dim, N);
    for i = 1 : N
        Ls(:, :, i) = chol(As(:, :, i), 'lower');
    end
    
    
    % Initial iterate for alternate descent method
    initalpha = rand(num, N); % ones(num, N) / num;
    Xinitial = reshape(Xalpha * max(pinv(initalpha), 0), dim, dim, num);
    
    % compute initial iterate
    fprintf('compute initial iterate\n');
    for i = 1 : 10
        pinvBs = pinv(reshape(Xinitial, dim * dim, []));
        initalpha = zeros(num, N);
        for j = 1 : N
%             a = norm(pinvBs * reshape(As(:, :, j), [], 1));
%             initalpha(:, j) = ones(num, 1) / sqrt(num) * a;
            initalpha(:, j) = max(pinvBs * reshape(As(:, :, j), [], 1), 0);
        end
%         i
%         initalpha
        Xinitial = reshape(Xalpha * max(pinv(initalpha), 0), dim, dim, num);
        
        fprintf('i:%d, Berr:%e, SCerr:%e\n', i, norm(Xinitial(:) - Xtrue(:)) / norm(Xtrue(:)), ...
            norm(initalpha(:) - alpha(:)) / norm(alpha(:)));
    end
    
    return;
% start the alternate descent method
%     alpha_i = alpha;
%     Dic_i = Xtrue;
    alpha_i = initalpha;
    Dic_i = Xinitial;
    maxiter = 50;
    lambdaB = logspace(-2, -6, maxiter);
    lambdaR = logspace(-2, -6, maxiter);

    for i = 1 : maxiter
        
%         Dic_i = reshape(Xalpha * max(pinv(alpha_i), 0), dim, dim, num); %%---
        % update the dictionary
        
%         fprintf('Dictionary learning subproblem\n');
        SolverParams.method = 'LRBFGS';
        SolverParams.Max_Iteration = 10;
        SolverParams.OutputGap = 100;
        SolverParams.LengthSY = 4;
        SolverParams.Num_pre_funs = 0;
        SolverParams.InitSteptype = 3;
        SolverParams.DEBUG = 0;
        SolverParams.Stop_Criterion = 2;
        HasHHR = 0;
        [XoptLRBFGS, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime] = TestSPDTensorDL(Ls, alpha_i, Dic_i, lambdaB(i), HasHHR, SolverParams);
        newDic_i = reshape(XoptLRBFGS.main, size(Dic_i));
        
        ftotal1 = f + sum(sum(alpha_i)) * lambdaR(i);

        
        
%         pinvBs = pinv(reshape(newDic_i, dim * dim, []));%%----
%         alpha_i = zeros(num, N);
%         for j = 1 : N
%             a = norm(pinvBs * reshape(As(:, :, j), [], 1));
%             alpha_i(:, j) = ones(num, 1) / sqrt(num) * a;
%         end%%----

%         fprintf('Sparse coding subproblem\n');
        SolverParams.method = 'LRBFGS';
        SolverParams.Max_Iteration = 10;
        SolverParams.OutputGap = 100;
        SolverParams.LengthSY = 4;
        SolverParams.Num_pre_funs = 0;
        SolverParams.InitSteptype = 3;
        SolverParams.DEBUG = 0;
        SolverParams.Stop_Criterion = 2;
        HasHHR = 0;
        [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times] = TestEucPosSpCd(Ls, newDic_i, alpha_i, lambdaR(i), HasHHR, SolverParams);
        newalpha_i = reshape(Xopt.main, size(alpha_i));
        
        ftotal2 = f;
        for j = 1 : num
            ftotal2 = ftotal2 + lambdaB(i) * trace(newDic_i(:, :, j));
        end
        
        ftotal3 = f - sum(sum(newalpha_i)) * lambdaR(i);
        
        fprintf('i:%d, Berr:%e, SCerr:%e, f1:%e, f2:%e, f3:%e\n', i, norm(newDic_i(:) - Xtrue(:)) / norm(Xtrue(:)), ...
            norm(newalpha_i(:) - alpha(:)) / norm(alpha(:)), ftotal1, ftotal2, ftotal3);
        Dic_i = newDic_i;
        alpha_i = newalpha_i;
    end
end

