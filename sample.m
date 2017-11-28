% Matlab implementation of FEM-H1 methodology
%
% for theory see
% PETSc user meeting, Vienna, Austria; June 28-30, 2016
% http://www.mcs.anl.gov/petsc/meetings/2016/posters/pospisil.pdf (poster)
%
% Pospisil L., Gagliardini P., Sawyer W., Horenko I.: On a scalable nonparametric denoising of time series signals
%

clear all

addpath 'problem'
addpath 'kmeans_FEMH1'

% generate some funny data
[x,T,K_sol,gamma_sol, mu_sol ] = problem1(0.3);
%[x,T,K_sol,gamma_sol, mu_sol ] = problem2(0.1);

% throw away K_sol, gamma_sol, mu_sol and try to find them

%% --- clustering with FEM-H1 regularisation ---
K1 = 3; % set number of clusters (should be computed for more values and right one chosen using AIC)
myeps = 10; % penalty-regularization parameter (should be computed for more values and right one chosen using L-curve)
[ mu2, gamma2, it2 ] = compute_kmeansh1(x,K1,myeps);

% compute reconstructed data
x2 = zeros(size(x))';
for k=1:K1
    x2 = x2 + gamma2((k-1)*T+1:k*T)*mu2(k);
end

% plot data
figure
subplot(K1,3,1:3:3*K1-2)
hold on
title('given data')
plot(1:T,x,'b','LineWidth',1.0)
xlabel('$t$', 'Interpreter', 'latex','FontSize',12)
ylabel('$x(t)$', 'Interpreter', 'latex','FontSize',12)
axis([0,T,min(x)-0.5,max(x)+0.5])
hold off

subplot(K1,3,3:3:3*K1)
hold on
title('reconstructed signal')
plot(1:T,x2,'Color',[0,0.5,0],'LineWidth',2.0)
xlabel('$t$', 'Interpreter', 'latex','FontSize',12)
ylabel('$x_{\mathrm{recovered}}(t)$', 'Interpreter', 'latex','FontSize',12)
axis([0,T,min(x)-0.5,max(x)+0.5])
hold off

% plot cluster indicator functions (gamma)
for i= 1:K1
    subplot(K1,3,2+3*(i-1))
    hold on
    plot(1:T,gamma2(((i-1)*T+1):i*T),'r','LineWidth',2.0)
    if i==K1
        xlabel('$t$', 'Interpreter', 'latex','FontSize',12)
    end
    ylabel(['$\gamma_{' num2str(i) '}(t)$'], 'Interpreter', 'latex','FontSize',12)
    axis([0,T,-0.2,1.2])
    hold off
end
