function [mu, x, xC, xN, xP] = find_best_fungi...
    (Ci, Ni, Pi, Ce, Ne, Pe, kappa, tau, delta, L, beta, sres, xres)

% sres is the number of iterative refinements of a solution
% that we run through for each value of x, and xres is the 
% number of different values for x that we try

% CC is the mass of carbon needed per unit volume of cell, given that
% maximal bias in the uptake of N instead of C
CC = max([Ci, Ni*Ce*delta/Ne]);
    
% hydrolases density cannot be so high that there is not enough local
% resource to synthesise them, which limits maximal rate of resource use,
% and relative rate of resource use must be less than L
U_max = min([L, (-1+(kappa^2)*Ce/CC)/tau, (-1+(kappa^2)*Ne/Ni)/tau]);

U_vector = linspace(U_max/xres, U_max, xres);
         
if (U_max <= 0) || ((kappa^2)*Pe/Pi < 1)
    
    mu = 0;
    x = 0;
    xC = 0; 
    xN = 0;
    xP = 0;
else
    
    mu_vector = zeros(xres, 1);
    xC_vector = zeros(xres, 1);
    xN_vector = zeros(xres, 1);
    xP_vector = zeros(xres, 1);
    
    U = U_vector(1);
    
    [mu_vector(1), xC_vector(1), xN_vector(1), xP_vector(1)] = ...
        find_mu_given_U_fungi(U, 0, 0, 0, 0, Ci, Ni, Pi, ...
            Ce, Ne, Pe, kappa, tau, delta, beta, xres);

    for i = 2:xres
        
        mu1 = mu_vector(i-1);
        xC1 = xC_vector(i-1);
        xN1 = xN_vector(i-1);
        xP1 = xP_vector(i-1);
        U = U_vector(i-1);
        
        [mu_vector(i), xC_vector(i), xN_vector(i), xP_vector(i)] = ...
            find_mu_given_U_fungi(U, mu1, xC1, xN1, xP1, Ci, Ni, Pi, ...
                Ce, Ne, Pe, kappa, tau, delta, beta, xres);
    end
    
    k = find(mu_vector == max(mu_vector), 1);
    
    if k == 1
        
        U_vector = linspace(U_vector(1)/xres, U_vector(2), xres);
        mu1 = mu_vector(1);
        xC1 = xC_vector(1);
        xN1 = xN_vector(1);
        xP1 = xP_vector(1);
        
    elseif k == xres
        
        U_vector = linspace(U_vector(xres-1), U_vector(xres), xres);
        mu1 = mu_vector(xres-1);
        xC1 = xC_vector(xres-1);
        xN1 = xN_vector(xres-1);
        xP1 = xP_vector(xres-1);
    else
        
        U_vector = linspace(U_vector(k-1), U_vector(k+1), xres);
        mu1 = mu_vector(k-1);
        xC1 = xC_vector(k-1);
        xN1 = xN_vector(k-1);
        xP1 = xP_vector(k-1);
    end
    
    for i = 1:xres
        
        [mu_vector(i), xC_vector(i),  xN_vector(i),  xP_vector(i)] = ...
            find_mu_given_U_fungi(U_vector(i), mu1, xC1, xN1, xP1,...
            Ci, Ni, Pi, Ce, Ne, Pe, kappa, tau, delta, beta, sres);
        
        mu1 = mu_vector(i);
        xC1 = xC_vector(i);
        xN1 = xN_vector(i);
        xP1 = xP_vector(i);
    end
    
    k = find(mu_vector == max(mu_vector), 1);
    
    mu = mu_vector(k);
    xC = xC_vector(k);
    xN = xN_vector(k);
    xP = xP_vector(k); 
    x = xC + xN + xP;
end
        