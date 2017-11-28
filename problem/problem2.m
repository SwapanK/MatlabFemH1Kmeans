function [ x, T, K, gamma, mu ] = problem2(sigma)
%PROBLEM1 set some funny values and generate problem
%
% x        generated data
% T        length of time-series
% K        number of clusters
% gamma    model indicator functions (0=inactive, 1=active)
% mu       mean values for each cluster
% sigma    noise parameter (standard deviation)
%
% Created by lukas.pospisil@usi.ch, Lugano, 2016

% length
T = 500;

% number of clusters
K = 3;

% mean values for clusters
mu{1} = -10;
mu{2} = 2;
mu{3} = 35;

% generate some funny jumping function
gamma{1} = zeros(1,T);
gamma{2} = zeros(1,T);
gamma{3} = zeros(1,T);

coeff = T/7;
for t = 1:T
   if t<=coeff
       gamma{1}(t) = 1;
   end
   if t>coeff && t <= 2*coeff
       gamma{2}(t) = 1;
   end
   if t>2*coeff && t <= 3*coeff
       gamma{3}(t) = 1;
   end
   if t>3*coeff && t <= 4*coeff
       gamma{1}(t) = 1;
   end  
   if t>4*coeff && t <= 5*coeff
       gamma{2}(t) = 1;
   end
   if t>5*coeff && t <= 6*coeff
       gamma{1}(t) = 1;
   end  
   if t>6*coeff
       gamma{2}(t) = 1;
   end  
end

% put data inside general generator
x = generate_kmeans(gamma,mu,sigma);

end

