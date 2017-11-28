function [ theta, gamma, it, Lit ] = compute_kmeansh1( x,K,myeps )
%COMPUTE_KMEANSH1 
%
% x        data
% K        number of clusters
% myeps    penalty-regularisation parameter
% theta    model parameters on each cluster
% gamma    model indicator functions
% it       number of iterations
%
% Created by lukas.pospisil@usi.ch, Lugano, 2016

disp('KMEANS FEM-H1 regularisation:')

% in this implementation I assume 1D data
T = length(x);

% generate initial random feasible gamma 
gamma = rand(K*T,1);
gamma_sum = zeros(T,1);
for k=1:K
    gamma_sum = gamma_sum+gamma((k-1)*T+1:k*T);
end
for k=1:K
    gamma((k-1)*T+1:k*T) = gamma((k-1)*T+1:k*T)./gamma_sum;
end

% create H (constant in outer it)
Hblock = 2*diag(ones(T,1)) - diag(ones(T-1,1),1) - diag(ones(T-1,1),-1);
Hblock(1,1) = 1;
Hblock(end,end) = 1;
Hblock = myeps*sparse(Hblock);
H = zeros(K*T,K*T);
for k = 1:K
   H((k-1)*T+1:k*T,(k-1)*T+1:k*T) = Hblock; 
end

% create equality constraints (constant in outer it)
B = zeros(T,K*T);
for k = 1:K
   B(:,(k-1)*T+1:k*T) = eye(T);
end
B = sparse(B);
c = ones(T,1);

% lower bound
l = zeros(K*T,1);

% this will be linear term in QP (changing in every outer it)
b = zeros(size(gamma));

% settings of algorithm (quadprog)
%options = optimoptions('quadprog');
options.Algorithm = 'interior-point-convex';
options.Display = 'none';

% here will be stored slution of model parameters - one for each cluster
theta = zeros(K,1);

% initial object function value
L = realmax;

it = 0; % iteration counter
while it < 1000 % practical stopping criteria after computing new L (see "break")

    % compute Theta
    for k=1:K
       if sum(gamma((k-1)*T+1:k*T)) ~= 0 % maybe gamma_k = 0 ?
           theta(k) = dot(gamma((k-1)*T+1:k*T),x)/sum(gamma((k-1)*T+1:k*T)); 
       else
           theta(k) = 0;
       end
    end

    % compute new gamma
    % prepare new linear term in QP
    for t = 1:T
        for k = 1:K
            b((k-1)*T+t) = (x(t) - theta(k))^2; % local model error
        end
    end
    % I would like to use gamma from previous outer iteration as initial
    % approximation in quadprog, but unfortunatelly ..
    [gamma, Lnew] = quadprog(H,b,[],[],B,c,l,[],gamma,options);

    % compute function value
    Lold = L;
    L = Lnew;

    disp([' it=' num2str(it) ', L=' num2str(L)]);
    
    if abs(L - Lold) < 10^-4
        break; % stop outer "while" cycle
    end
    
    it = it + 1;
    
    Lit(it) = L; % for postprocessing
end    


end

