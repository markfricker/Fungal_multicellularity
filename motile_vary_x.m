% xres is the number of different relative densities of 
% exoenzymes that we try
xres = 1000;

% sres is the number of iterative refinements of a solution
% that we run through for each value of x, 
sres = 50;

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
Ci = 0.165;
Ni = 0.032;
Pi = 0.005;

% the carbon use efficiency is the fraction of C for growth versus total C
% for respiration and growth
CUE = 0.5;

% Ct is the total C required, including respiratory C set by the carbon use
% efficiency
Ct = Ci/CUE;

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
% relative to Ct + Ni + Pi
L = lambda/(Ct + Ni + Pi);

M_tot = 12*C_to_P + 14*N_to_P + 31;
Ce = density*12*C_to_P/M_tot;
Ne = density*14*N_to_P/M_tot;
Pe = density*31/M_tot;

% phi is the correction term
phi = (Ct + Ni)/(Ct + Ni + Pi);

% CC is the mass of carbon needed per unit volume of cell, given that
% maximal bias in the uptake of N instead of C
CC = max([Ct, Ni*Ce*delta/Ne]);
    
% theta is the fraction of digested resource that is strictly  
% necessary for growth, when everything is used for hydrolases.
% this is the lower bound, which gives maximal x
theta = (Ni + Ct)/(Ni + CC);

% optimal density of hydrolases x <= L*tau/(phi*theta), 
% and mass of hydrolases cannot be greater than mass of locally 
% present resource, so only consider relevant range
x_max = min([L*tau/(phi*theta), (kappa^2)*Ce/CC]);

x_vector = linspace(x_max/xres, x_max, xres);
M_vector = zeros(xres, 1);

for i = 1:xres
    
    x = x_vector(i);
    
    % theta is the fraction of digested resource that is strictly  
    % necessary for growth
    theta = (Pi + (Ni + Ct)*(1 + alpha + x))/...
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
    mu = eta;
    
    for dum = 1:sres
        
        mu = eta*(exp(mu*omega/eta) - 1)/...
            (exp(mu*omega/eta) - exp(-mu*tau/theta));
    end
    
    if (isreal(mu) == 1) && (mu > 0)
        
        M_vector(i) = mu;
    end
end

[M, x] = find_best_motile...
    (Ct, Ni, Pi, Ce, Ne, Pe, kappa, tau, delta, L, alpha, sres, xres);

figure(1)
plot(x_vector, M_vector)
hold on
plot(x, M, 'o')
xlabel('Relative density of exoenzymes')
ylabel('Apparent growth rate')