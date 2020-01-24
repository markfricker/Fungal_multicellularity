function [mu, eta, x] = find_best_autolytic...
    (Ct, Ni, Pi, Ce, Ne, Pe, kappa, tau, delta, L, epsilon, sres, xres)

% sres is the number of iterative refinements of a solution
% that we run through for each value of x, and xres is the 
% number of different values for x that we try

% phi is the correction term
phi = (Ct + Ni)/(Ct + Ni + Pi);

% CC is the mass of carbon needed per unit volume of cell, given that
% maximal bias in the uptake of N instead of C
CC = max([Ct, Ni*Ce*delta/Ne]);
    
% theta is the fraction of digested resource that is strictly  
% necessary for growth, when everything is used for hydrolases
theta = (Ni + Ct)/(Ni + CC);
    
% optimal density of hydrolases x <= L*tau/(phi*theta), 
% and mass of hydrolases cannot be greater than mass of locally 
% present resource, as at best recycled resource pays for the
% cell itself, so only consider relevant range     
x_max = min([L*tau/(phi*theta), ...
    (kappa^2)*Ce/CC, (kappa^2)*Ne/Ni]);
      
x_vector = linspace(x_max/xres, x_max, xres);

if (x_max <= 0)
    
    mu = 0;
    eta = 0;
    x = 0;
else
    
    mu_vector = zeros(xres, 1);
    eta_vector = zeros(xres, 1);
    x = x_vector(1);
    
    [mu_vector(1), eta_vector(1)] = find_mu_given_x_autolytic(x, 0, 0, ...
        Ct, Ni, Pi, Ce, Ne, Pe, kappa, tau, L, delta, epsilon, sres);
    
    for i = 2:xres
        
        x = x_vector(i);
        mu1 = mu_vector(i-1);
        eta1 = eta_vector(i-1);
        
        [mu_vector(i), eta_vector(i)] = find_mu_given_x_autolytic...
            (x, mu1, eta1, Ct, Ni, Pi, Ce, Ne, Pe, ...
            kappa, tau, L, delta, epsilon, sres);
        
    end
    
    k = find(mu_vector == max(mu_vector), 1);
    
    if k == 1
        
        x_vector = linspace(x_vector(1)/xres, x_vector(2), xres);
        mu1 = mu_vector(1);
        eta1 = eta_vector(1);
        
    elseif k == xres
        
        x_vector = linspace(x_vector(xres-1), x_vector(xres), xres);
        mu1 = mu_vector(xres-1);
        eta1 = eta_vector(xres-1);
    else
        x_vector = linspace(x_vector(k-1), x_vector(k+1), xres);
        mu1 = mu_vector(k-1);
        eta1 = eta_vector(k-1);
    end
    
    mu_vector = zeros(xres, 1);
    eta_vector = zeros(xres, 1);
    
    for i = 1:xres
        
        x = x_vector(i);
        
        [mu_vector(i), eta_vector(i)] = find_mu_given_x_autolytic...
            (x, mu1, eta1, Ct, Ni, Pi, Ce, Ne, Pe, ...
            kappa, tau, L, delta, epsilon, sres);

        mu1 = mu_vector(i);
        eta1 = eta_vector(i);
    end
    
    k = find(mu_vector == max(mu_vector), 1);
    mu = mu_vector(k);
    eta = eta_vector(k);
    x = x_vector(k);  
end
        