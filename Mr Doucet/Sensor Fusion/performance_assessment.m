clear; close all; clc;

nb_tests = 50;
filter = 'PF';
sigma_P2 = 12000;
randomness = false;
anisotropic = true;

indices = zeros(4,nb_tests);
parfor i=1:nb_tests
    disp(i)
    [CQ,AI,covX,covY] = runFilter(filter,sigma_P2,randomness,anisotropic);
    indices(:,i)=[CQ;AI;covX;covY];
end

% figure
% histogram(indices(1,:),100);
% title('Repartition of CQ over 500 experiments','interpreter','latex')

figure
[fCQ,CQ] = ksdensity(indices(1,:));
plot(CQ,(nb_tests/sum(fCQ))*fCQ);
title('Repartition of CQ over 50 experiments','interpreter','latex')
% if anisotropic
%     printEPS(['CQ_' filter '_aniso']);
% else
%     printEPS(['CQ_' filter '_iso']);
% end

% figure
% histogram(indices(2,:),100);
% title('Repartition of AI over 500 experiments','interpreter','latex')

figure
[fAI,AI] = ksdensity(indices(2,:));
plot(AI,(nb_tests/sum(fAI))*fAI);
title('Repartition of AI over 50 experiments','interpreter','latex')
% if anisotropic
%     printEPS(['AI_' filter '_aniso']);
% else
%     printEPS(['AI_' filter '_iso']);
% end
[m,i] = max(fAI);
disp(['AIrep' AI(i)])

figure
subplot(2,1,1)
[fx,x] = ksdensity(indices(3,:));
plot(x,(nb_tests/sum(fx))*fx);
title('Repartition of final covariance over 50 experiments','interpreter','latex')
xlabel('cov(x)','interpreter','latex')
subplot(2,1,2)
[fy,y] = ksdensity(indices(3,:));
plot(y,(nb_tests/sum(fy))*fy);
xlabel('cov(y)','interpreter','latex')
% if anisotropic
%     printEPS(['cov_' filter '_aniso']);
% else
%     printEPS(['cov_' filter '_iso']);
% end