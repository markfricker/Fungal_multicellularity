function [mu] = find_mu_given_x_immobile...
    (x, mu1, Ci, Ni, Pi, Ce, Ne, Pe, kappa, tau, L, delta, res)

% mu1 is the initial guess for what the growth rate will be

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

% eta is the specific growth rate of immobile cells
eta = min([L/(1 + phi*x), (theta*phi*x/tau)/(1 + phi*x)]);

% mu is our initial guess for the growth rate of persistent
% cells when the rate of digestion is D
if mu1 > 0
    
    mu = mu1;
else
    mu = eta;
end

for dum = 1:res
    
    mu = eta*(1 - exp(-mu*omega/eta));
end

if (isreal(mu) ~= 1) || (mu < 0) || (isnan(mu) == 1)
    
    mu = 0;
end
