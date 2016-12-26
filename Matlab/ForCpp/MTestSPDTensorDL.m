function MTestSPDTensorDL()
    r = floor(rand() * 10000000)
%     r = 268764
%     r = 7740966
    rand('state', r);
    randn('state', r);
    dim = 4;%% 20;
    num = 10;%% 100;
    N = 200;%% 500;
    numatoms = 5;
    
%     ftrue = 0;
    Xtrue = zeros(dim, dim, num);
    for i = 1 : num
        tmp = randn(dim, dim * 10);
        Xtrue(:, :, i) = tmp * tmp';
%         ftrue = ftrue  + trace(Xtrue(:, :, i));
    end
%     ftrue
    
    alpha = zeros(num, N);
    for i = 1 : N
        values = rand(numatoms, 1);
        idx = randperm(num, numatoms);
        for j = 1 : numatoms
            alpha(idx(j), i) = values(j);
        end
    end
%     alpha
%     svd(alpha)
%     sum(alpha')
    Xalpha = reshape(Xtrue, dim * dim, num) * alpha;
    As = reshape(Xalpha, dim, dim, N);
    
    Ls = zeros(dim, dim, N);
    for i = 1 : N
        Ls(:, :, i) = chol(As(:, :, i), 'lower');
    end
    
    Xinitial = reshape(Xalpha * max(pinv(alpha), 0), dim, dim, num);
    
%     % initial iterate
%     Xinitial = zeros(dim, dim, num);
%     for i = 1 : num
%         tmp = randn(dim, dim);
%         Xinitial(:, :, i) = tmp * tmp'* 1e-0;% + Xtrue(:, :, i);
% %         Xinitial(:, :, i) = eye(dim);
%     end
    lambdaB = 0;

    SolverParams.method = 'LRBFGS';
%     SolverParams.method = 'RTRSR1';
%     SolverParams.method = 'RTRNewton';
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
    fprintf('start LRBFGS\n');
    [XoptLRBFGS, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime] = TestSPDTensorDL(Ls, alpha, Xinitial, lambdaB, HasHHR, SolverParams);

    
%     Xsolu = reshape(Xopt.main, size(Xtrue));
%     trX = 0;
%     for i = 1 : size(Xsolu, 3)
%         trX = trX + sum(diag(Xsolu(:, :, i)));
%     end
%     trX
    
%     Xtrue
%     
%     norm(Xopt.main - reshape(Xtrue, [], 1)) / norm(reshape(Xtrue, [], 1))
    
    SolverParams.method = 'RCG';
    SolverParams.IsCheckParams = 1;
    SolverParams.Max_Iteration = 500;
    SolverParams.OutputGap = 100;
    SolverParams.InitSteptype = 1;
    SolverParams.Stop_Criterion = 2;
    SolverParams.DEBUG = 2;
    SolverParams.RCGmethod = 1;
    HasHHR = 0;
    fprintf('start RCG\n');
    [XoptRCG, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime] = TestSPDTensorDL(Ls, alpha, Xinitial, lambdaB, HasHHR, SolverParams);

    norm(XoptLRBFGS.main - XoptRCG.main) / norm(XoptLRBFGS.main)
%     SolverParams.method = 'LRTRSR1';
%     SolverParams.IsCheckParams = 1;
%     SolverParams.Max_Iteration = 500;
%     SolverParams.OutputGap = 50;
%     SolverParams.InitSteptype = 1;
%     SolverParams.Stop_Criterion = 1;
%     SolverParams.DEBUG = 2;
%     HasHHR = 0;
%     [Xopt, f, gf, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime] = TestSPDTensorDL(Ls, alpha, Xinitial, 1, HasHHR, SolverParams);
end
% 
% %% generate test data set
% function As = matrix_data(n, k, data_type)
%     As = cell(k,1);
%     if(data_type == 1) %ill
%         CN = 4;
%         for i = 1 : CN
%             [O,~] = qr(randn(n));
%             f = 0;
%             D = diag([[rand(1,n-1)+1],10^(-f)]);
%             As{i} = O * D * O';
%         end
%         for i = (CN+1) : k
%             [O,~] = qr(randn(n));
%             f = 5;
%             %D = 100000 * diag([[rand(1,n-1)+1],10^(-f)]);
%             D = diag([[rand(1,n-1)+1],10^(-f)]);
%             As{i} = O * D * O';
%             % aha = cond(As{i})
%         end
%     end
% 
% 
%     if(data_type == 3)
%         [O,~] = qr(randn(n));
%         f = 0;
%         D = diag([[rand(1,n-1)+1],10^(-f)]);
%         As{1} = O * D * O';
% 
%         for i = 2 : k
%             [O,~] = qr(randn(n));
%             f = 5;
%             m = 70;
%             Dx = diag([rand(1,m)+1, (rand(1,n-m)+1)*10^(f)]);
%             ill = O * Dx * O';
%             As{i} = As{1} * ill * As{1};
%             % aha = cond(As{i})
%         end
% 
%     end
% 
% 
% 
%     if(data_type == 2)%structured well
%         grp = floor(k/3);
%         for i = 1 : grp
%             f = 0;
%             [O,~] = qr(randn(n));
%             D = diag([[rand(1,n-1)+1],10^(-f)]);
%             As{i} = O * D * O';
%         end
%         for i = (grp+1) : (2*grp)
%             f = 1;
%             [O,~] = qr(randn(n));
%             D = diag([[rand(1,n-1)+1],10^(-f)]);
%             As{i} = O * D * O';
%         end
%         for i = (2*grp + 1) : k
%             f = 2;
%             [O,~] = qr(randn(n));
%             D = diag([[rand(1,n-1)+1],10^(-f)]);
%             As{i} = O * D * O';
%         end
%     end
% end
