function [ x ] = generate_kmeans(gamma,mu,sigma)
%GENERATE_KMEANS generate testing data for non-stationary mean problem
%
% x        generated data
% gamma    model indicator functions (0=inactive, 1=active)
% mu       mean values for each cluster
% sigma    noise parameter (standard deviation)
%
% Created by lukas.pospisil@usi.ch, Lugano, 2016

% length of time-series
T = length(gamma{1});

% dimension of data
xdim = length(mu{1});

% number of clusters
K = length(gamma);

x = zeros(xdim,T);
for t = 1:T
    % get id of cluster from provided gamma
    this_k = 0;
    for k=1:K
        if gamma{k}(t) == 1
            this_k = k;
        end
    end
    
    % constant term
    x(:,t) = mu{this_k};

    % add noise
    x(:,t) = x(:,t) + sigma*randn(xdim,1);
end

end

