function [mu] = find_mu_given_x_motile...
    (x, mu1, Ci, Ni, Pi, Ce, Ne, Pe, kappa, tau, delta, alpha, res)

% mu1 is the initial guess for what the growth rate will be

% phi is the correction term
phi = (Ci + Ni)/(Ci + Ni + Pi);

% CC is the mass of carbon needed per unit volume of cell, given that
% maximal bias in the uptake of N instead of C
CC = max([Ci, Ni*Ce*delta/Ne]);
    
% theta is the fraction of digested resource that is strictly  
% necessary for growth
theta = (Pi + (Ni + Ci)*(1 + alpha + x))/...
        (Pi + (Ni + CC)*(1 + alpha + x));
    
% Omega is the number of duplicates that can be made of the cell and
% its hydrolases, given the local supply of nutrients
omega = min([((kappa^2)*Ce - CC*x)/(CC*(1 + alpha + x)), ...
              ((kappa^2)*Ne - Ni*x)/(Ni*(1 + alpha + x)), ...
              (kappa^2)*Pe/Pi]);
    
% eta is the specific growth rate of motile cells
eta = (theta*phi*x/tau)/(1 + phi*alpha + phi*x);
    
% mu is our initial guess for the growth rate of mobile
% cells when the relative density of exoenzymes is x 
if mu1 == 0 
    
    mu = eta;
else
    mu = mu1;
end
        
for dum = 1:res
        
    mu = eta*(exp(mu*omega/eta) - 1)/...
            (exp(mu*omega/eta) - exp(-mu*tau/theta));
end
    
if (isreal(mu) ~= 1) || (mu < 0) || (isnan(mu) == 1)
        
    mu = 0;
end