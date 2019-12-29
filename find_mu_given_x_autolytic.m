function [mu, eta] = find_mu_given_x_autolytic(x, mu1, eta1, ...
    Ci, Ni, Pi, Ce, Ne, Pe, kappa, tau, L, delta, epsilon, res)

% mu1 is the initial guess for what the growth rate will be
% eta1 is the initial guess for the specific growth rate

% phi is the correction term
phi = (Ci + Ni)/(Ci + Ni + Pi);

% CC is the mass of carbon needed per unit volume of cell, given that
% maximal bias in the uptake of N instead of C
CC = max([Ci, Ni*Ce*delta/Ne]);

% theta is the fraction of digested resource that is strictly
% necessary for growth
theta = (Pi + (Ni + Ci)*(1 + x))/(Pi + (Ni + CC)*(1 + x));

% omega is the number of duplicates that can be made of the cell and
% its hydrolases, given the local supply of nutrients
omega = min([(kappa^2)*Ce/(CC*(1 + x)), ...
             (kappa^2)*Ne/(Ni*(1 + x)), (kappa^2)*Pe/Pi]);

% mu is our initial guess for the growth rate of autolytic cells
if mu1 == 0
    
    mu = (theta*phi*x/tau);
else
    mu = mu1;
end

% eta is our initial guess for the specific growth rate
% of autolytic cells
if eta1 == 0
    
    eta = mu;
else
    eta = eta1;
end

for dum = 1:res
    
    mu = eta*(1 - exp(-mu*omega/eta));
    
    U = (theta*phi*x/tau) + epsilon*mu/(exp(mu*omega/eta) - 1);
    
    eta = min(U, L)/(1 + phi*x);
end

if (isreal(mu) ~= 1) || (mu < 0) || ...
    (isnan(mu) == 1) || exp(-mu*omega/eta) > 0.99
    
    mu = 0;
    eta = 0;
end