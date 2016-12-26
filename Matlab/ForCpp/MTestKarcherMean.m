function MTestKarcherMean(C)

    if(nargin == 0)
        Load_Datasets;
        C = CM;
    end
    
    fprintf('Test for Karcher Mean with Path Straightening in Shape Space. \n');
    numP = 101;
    dim = 2;
    numS = 4;
    numC = 11;
    ncluster = 1;
    
    for i = 1 : numS
        Cs(:, :, i) = [C(:, :, i+20*ncluster), C(:, 1, i+20*ncluster)];
    end
        
    Cinitial = [C(:, :, 80 + 20 * ncluster), C(:, 1, 21 + 20 * ncluster)];


%     for i = 2 : numS
%         tempC = updateC(Cs(:,:,1), Cs(:,:,i));
%         Cs(:,:,i) = tempC;
%         figure(i);clf
%         size(Cs(:,:,1), 2)
%         plot_matching_curves(Cs(:,:,1), Cs(:,:,i), linspace(0, 1, size(Cs(:,:,1), 2)), 1, '');
%         hold on
%     end

     for i  = 2 : numS
        tempC = updateC(Cs(:,:,1), Cs(:,:,i));
        Cs(:,:,i) = tempC;
        figure(i);clf
        size(Cs(:,:,i), 2)
        plot_matching_curves(Cs(:,:,1), Cs(:,:,i), linspace(0, 1, size(Cs(:,:,1), 2)), 1, '');
        hold on 
     end 
    
% %     %%=================== Get Initial =============================
% % %     rotated = 1;
% % %     isclosed = 1;
% % %     onlyDP = 0;
% % %     auto = 2;
% % %     w = 0.01;
% % %     for i = 1 : numS
% % % 
% % %         C1 = Cs(:,:,i);
% % %         C2 = Cinitial;
% % %         clear mex
% % %         [XQNopt, swap, QNfopts, QNcomtime] = TestElasticCurvesRO(C1, C2, w, rotated, isclosed, onlyDP, 4, 'LRBFGS', auto);
% % %         QNm = XQNopt(end);
% % %         O = XQNopt(end - 4 : end - 1);
% % %         QNO = reshape(0, 2, 2);
% % %         QNgamma = XQNopt(1 : end - 5);
% % %         if (swap)
% % %             QNgamma = invertGamma(QNgamma);
% % %             QNO = QNO';
% % %             QNm = 1 - QNm;
% % %         end 
% % %         gammaQN = mod(QNgamma + QNm, 1);
% % %         plot_matching_curves(C1', (C2 * QNO')', gammaQN, 1, strcat('R, L:', num2str(QNfopts(1))));
% % % 
% % %      end 
% %     
    figure(100);clf
    for i = 1 : numS
        subplot(1, numS, i);
        plot(Cs(1,:,i), Cs(2,:,i), '.');
        axis equal
        axis off
    end
    pause(0.1)

    for i = 1 : numS
        qs(:,:,i) = curve_to_q(Cs(:,:,i))';
    end 
 

    q_initial = curve_to_q(Cinitial)'; %qs{1};
    %q_initial
    %size(q_initial)
    figure(101);clf
    plot(Cinitial(1, :), Cinitial(2, :), '.');
    axis equal
    axis off
    pause(0.1);
     
     [Mean] = TestKarcherMean(qs, q_initial, numP, dim, numS, numC);
    CMean = q_to_curve(Mean');   
    figure(200);
    plot(CMean(1,:), CMean(2,:), '.');
    axis equal
    axis off

% Mean = q_to_curve(Cmean');
%     figure(200);
%     plot(Mean(1,:), Mean(2,:), '.');
%     axis equal
%     axis off

end 


function C2 = updateC(C1, C2)
    N = size(C1, 2);

    rotated = 1;
    isclosed = 1;
    auto = 2;
    onlyDP = 0;
    clear mex
    [XQNopt, swap2, QNfopts, QNcomtime] = TestElasticCurvesRO(C1', C2', 0.01, rotated, isclosed, onlyDP, 4, 'LRBFGS', auto);
    QNm = XQNopt(end)
    O = XQNopt(end - 4 : end - 1);
    QNO = reshape(O, 2, 2);
    QNgamma = XQNopt(1 : end - 5);
    if(swap2)
        QNgamma = invertGamma(QNgamma);
        QNO = QNO';
        QNm1 = QNm %%---
        QNm = 1 - QNm;
        QNm2 = QNm %%---
    end
    %         figure(1);clf
    x1.m = QNm;
    x1.O = QNO;
    x1.l = reshape(sqrt(gradient(QNgamma) / N), [], 1);
    x1.l = x1.l / sqrt(trapz(linspace(0, 1, length(x1.l)), x1.l .* x1.l));
    %         subplot(1, 3, 1); hold on
    QNgamma = mod(QNgamma + QNm, 1);
    %         scatter(idx, QNgamma, '.', 'k');
    %         axis([0, 1, 0, 1]);

    %     subplot(1, 3, 2);
    %     plot_matching_curves(C1, (QNO * C2), QNgamma, 1, 'Initial');
    %     hold on
    %     pause(0.01);
    %     subplot(1, 3, 3);
    C2 = updateC2((QNO * C2), QNgamma);
    C2(:, end) = C2(:, 1);
end

function C2 = updateC2(C2, gamma)
    N = size(C2, 2);
    ppC2 = spline((0:(N-1))/(N-1), C2);
    for i = 1 : N
        VV = ppval(ppC2, gamma(i));
        C2(:, i) = [VV(1); VV(2)];
    end
end

function gamI = invertGamma(gam)
    N = length(gam);
    x = [0:(N-1)]/(N-1);
    gamI = interp1(gam,x,x,'linear');
end

function output = Center_C(C)
    [n, N] = size(C);
    ave = sum(C')' / N;
    output = C;
    for i = 1 : N
        output(:, i) = output(:, i) - ave;
    end
end

function plot_matching_curves(C1, C2, gamma, L, titlestr)
    C1 = Center_C(C1);
    C1 = C1 / norm(C1);
    C2 = Center_C(C2);
    C2 = C2 / norm(C2);
    C2(1, :) = C2(1, :) - min(C2(1, :)) + max(C1(1, :)) + 0.02;
    N = size(C1, 2);
    color_hsv = hsv(N);
    plot(C1(1, :), C1(2, :), 'b')
    hold on
    plot(C2(1, :), C2(2, :), 'r')
    title(titlestr)
    axis equal
    axis tight
    ppC2 = spline((0:(N-1))/(N-1) * L, C2);
    for i = 1 : N
        VV = ppval(ppC2, gamma(i));
        scatter([C1(1, i), VV(1)], [C1(2, i), VV(2)], 10, [color_hsv(i, :)]);
    %         plot([C1(1, i), VV(1)], [C1(2, i), VV(2)], 'k');
        if(i == 1)
            scatter([C1(1, i), VV(1)], [C1(2, i), VV(2)], '*', 'k');
        end
    end
    axis off
end

function [q] = curve_to_q(p)                %SRVF
    [n,N] = size(p);                        %n dimension, N number of points
    xslope = (p(1, 2) - p(1, end - 1)) / (2 / (N - 1));   %
    yslope = (p(2, 2) - p(2, end - 1)) / (2 / (N - 1));
    pp = zeros(n, N + 2);
    pp(1, 1) = xslope;
    pp(1, end) = xslope;
    pp(2, 1) = yslope;
    pp(2, end) = yslope;
    pp(:, 2 : end - 1) = p;
    
    ppp = spline((0 : (N - 1)) / (N - 1), pp);     %spline: matlab���������ֵ����(��һ��߽�����)
    ppdp = first_deri(ppp);
    v = ppval(ppdp, (0 : (N - 1)) / (N - 1));      %y=ppval(pp,xx)��ֵ�����ϣ�xx���Ӧ��yy
    for i = 1:N
        L(i) = sqrt(norm(v(:,i),'fro'));
        if L(i) > 0.000001
            q(:,i) = v(:,i)/L(i);
        else
            q(:,i) = 0*ones(n,1);
        end
    end
    q = q/sqrt(InnerProd_Q(q,q));
end

function [p] = q_to_curve(q)
    [n,T] = size(q);
    for i = 1:T
        qnorm(i) = norm(q(:,i),'fro');  %Frobenius norm
    end
    for i = 1:n
        p(i,:) = [ cumtrapz( q(i,:).*qnorm )/(T) ] ; 
    end
end

function output = first_deri(ppf)
    ppf.coefs = ppf.coefs(:, 1 : 3);
    ppf.coefs(:, 1) = ppf.coefs(:, 1) * 3;
    ppf.coefs(:, 2) = ppf.coefs(:, 2) * 2;
    ppf.order = 3;
    output = ppf;
end

function val = InnerProd_Q(q1,q2)  %trapezoidal rule
    [n,T] = size(q1);
    val = trapz(linspace(0,1,T),sum(q1.*q2)); %
end
