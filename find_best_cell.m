function [M, x] = find_best_cell(Ci, Ni, Pi, Ce, Ne, tau, delta, L)

% res1 is the number of different values for x we try
res = 1000;

% phi is the correction term
phi = (Ci + Ni)/(Ci + Ni + Pi);

% CC is the mass of carbon needed per unit volume of cell, given that
% maximal bias in the uptake of N instead of C
CC = max([Ci, Ni*Ce*delta/Ne]);

% theta0 is the fraction of digested resource that is strictly  
% necessary for growth, when everything there are no hydrolases
theta0 = (Ci + Ni + Pi)/(CC + Ni + Pi);

% theta1 is the fraction of digested resource that is strictly  
% necessary for growth, when everything is used for hydrolases
theta1 = (Ni + Ci)/(Ni + CC);

% growth is maximal when x is just large to supply the cell so it can grow
% at the maximal rate. This depends on theta, which depends on x
x_vector = linspace(L*tau/(phi*theta0), L*tau/(phi*theta1), res);
M_vector = zeros(res, 1);

for i = 1:res
    
    x = x_vector(i);
    
    % theta is the fraction of digested resource that is strictly  
    % necessary for growth
    theta = (Pi + (Ni + Ci)*(1 + x))/(Pi + (Ni + CC)*(1 + x));
    
    mu = (theta*phi*x/tau)/(1 + phi*x);
    
    M_vector(i) = mu;
end

k = find(M_vector == max(M_vector), 1);

M = M_vector(k);
x = x_vector(k);