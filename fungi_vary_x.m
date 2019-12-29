% xres is the number of different relative densities of 
% exoenzymes that we try
xres = 1000;

% sres is the number of iterative refinements of a solution
% that we run through for each value of x, 
sres = 100;

% C_to_N and C_to_P ratio of environment
C_to_P = 2000;
N_to_P = 10;

% dry weight in grams per ml of the substrate
density = 0.5;

% kappa is the relative radius of the resource acquisition zone
kappa = 1;

% tau is the time in hours for hydrolases to digest their own mass
tau = 20;

% delta is the minimal fraction of C that must be digested in order to 
% digest the available N
delta = 0;

% Ci is the mass of carbon needed per unit volume of growth, in g per ml
Ci = 0.33;
Ni = 0.032;
Pi = 0.005;

% epsilon is the efficiency of recycling for autolytic cells
epsilon = 0.5;

% alpha is the mass of machinery needed for cell mobility, relative
% to the mass of essential machinery
alpha = 0.02;

% beta is the mass of material in vesicles, relative to the mass of
% the rest of the fungus, including hydrolases
beta = 0.1;

% lambda is the maximum rate of resource use per unit volume,
% in g per ml per hour
lambda = 0.3;

% L is the maximum rate of resource use per unit volume,
% relative to Ci + Ni + Pi
L = lambda/(Ci + Ni + Pi);

M_tot = 12*C_to_P + 14*N_to_P + 31;
Ce = density*12*C_to_P/M_tot;
Ne = density*14*N_to_P/M_tot;
Pe = density*31/M_tot;

% phi is the correction term
phi = (Ci + Ni)/(Ci + Ni + Pi);

% CC is the mass of carbon needed per unit volume of cell, given that
% maximal bias in the uptake of N instead of C
CC = max([Ci, Ni*Ce*delta/Ne]);
    
% hydrolases density cannot be so high that there is not enough local
% resource to synthesise them, which limits maximal rate of resource use,
% and relative rate of resource use must be less than L
U_max = min([L, (-1+(kappa^2)*Ce/CC)/tau, (-1+(kappa^2)*Ne/Ni)/tau]);

U_vector = linspace(U_max/xres, U_max, xres);

M_vector = zeros(xres, 1);
xC_vector = zeros(xres, 1);
xN_vector = zeros(xres, 1);
xP_vector = zeros(xres, 1);
x_vector = zeros(xres, 1);

for i = 1:xres
    
    % find growth rate for given rate of synthesis, guessing the
    % growth rate and rates of digestion and updating iteratively
    U = U_vector(i);
    
    % xC is our initial guess for the relative density of exoenzymes, 
    % given the rate of synthesis S. Likewise for xN and xP
    xC = U*tau*CC/(Ci + Ni + Pi);
    xN = U*tau*Ni/(Ci + Ni + Pi);
    xP = U*tau*Pi/(Ci + Ni + Pi);
    
    x = xC + xN + xP;
    
    % mu is our initial guess for the growth rate of fungi
    mu = U/(1 + phi*(x + beta + beta*x));
    
    for dum = 1:sres
        
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

    if (isreal(mu) == 1) && (mu > 0) && (x <= x_max)
        
        M_vector(i) = mu;
        xC_vector(i) = xC;
        xN_vector(i) = xN;
        xP_vector(i) = xP;
        x_vector(i) = xC + xN + xP;
    else
        x_vector(i) = max(x_vector);
    end
end

[M, x, xC, xN, xP] = find_best_fungi...
     (Ci, Ni, Pi, Ce, Ne, Pe, kappa, tau, delta, L, beta, sres, xres);

figure(1)
plot(x_vector, M_vector)
hold on
plot(x,M,'o')
xlabel('Relative density of exoenzymes')
ylabel('Apparent growth rate')