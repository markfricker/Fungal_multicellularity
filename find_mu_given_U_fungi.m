function [mu, xC, xN, xP] = find_mu_given_U_fungi...
    (U, mu1, xC1, xN1, xP1, Ci, Ni, Pi, Ce, Ne, Pe, ...
    kappa, tau, delta, beta, res)

% mu1 is the initial guess for the apparent growth rate of the colony
% xC1 is the initial guess for the relative density of C digesting
% exoenzymes, and likewise for xN1 and xP1

% phi is the correction term
phi = (Ci + Ni)/(Ci + Ni + Pi);

% CC is the mass of carbon needed per unit volume of cell, given that
% maximal bias in the uptake of N instead of C
CC = max([Ci, Ni*Ce*delta/Ne]);

if xC1 == 0
    
    xC = U*tau*CC/(CC + Ni + Pi);
else
    xC = xC1;
end
if xN1 == 0
    
    xN = U*tau*Ni/(CC + Ni + Pi);
else
    xN = xN1;
end
if xP1 == 0
    
    xP = U*tau*Pi/(CC + Ni + Pi);
else
    xP = xP1;
end

x = xC + xN + xP;

if mu1 == 0
    
    % mu is our initial guess for the growth rate of fungi
    mu = U/(1 + phi*(x + beta + beta*x));
else
    mu = mu1;
end

for dum = 1:res
    
    % fC is the mass fraction of required resource that is C
    fC = Ci*(1 + x)*(1 + beta)/...
        (Pi + (Ci + Ni)*(1 + x)*(1 + beta));
    
    % fN is the mass fraction of required resource that is N
    fN = Ni*(1 + x)*(1 + beta)/...
        (Pi + (Ci + Ni)*(1 + x)*(1 + beta));
    
    % fP is the mass fraction of required resource that is P
    fP = Pi/(Pi + (Ci + Ni)*(1 + x)*(1 + beta));
    
    % tC is the time taken to exhaust the local supply of C,
    % and likewise for tN and tP
    tC = (kappa^2)*Ce*tau/(xC*(Ci + Ni + Pi));
    tN = (kappa^2)*Ne*tau/(xN*(Ci + Ni + Pi));
    tP = (kappa^2)*Pe*tau/(xP*(Ci + Ni + Pi));
    
    % xP is the relative density of P digesting exoenzymes
    xP = fP*U*tau/(1 - exp(-mu*tP));
    
    % xN is the relative density of N digesting exoenzymes
    xN = fN*U*tau/(1 - exp(-mu*tN));
    
    % xC is the relative density of C digesting exoenzymes
    xC = max([fC*U*tau/(1 - exp(-mu*tC)), xN*Ce*delta/Ne]);
    
    x = xC + xN + xP;
    
    mu = U/(1 + phi*(x + beta + beta*x));
end

 x_max = min([-1+(kappa^2)*Ce/(CC*(1+beta)), ...
              -1+(kappa^2)*Ne/(Ni*(1+beta))]);
             
if (isreal(mu) ~= 1) || (mu < 0) || (isnan(mu) == 1) || (x > x_max)
    
    mu = 0;
    xC = 0;
    xN = 0;
    xP = 0;
end
