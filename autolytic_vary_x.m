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
kappa = 5;

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
    
% theta is the fraction of digested resource that is strictly  
% necessary for growth, when everything is used for hydrolases
theta = (Ni + Ci)/(Ni + CC);
 
x_max = min([L*tau/(phi*theta), ...
    -1+(kappa^2)*Ce/CC, -1+(kappa^2)*Ne/Ni]);

x_vector = linspace(x_max/xres, x_max, xres);
M_vector = zeros(xres, 1);

for i = 1:xres
    
    x = x_vector(i);
    
    % theta is the fraction of digested resource that is strictly  
    % necessary for growth
    theta = (Pi + (Ni + Ci)*(1 + x))/(Pi + (Ni + CC)*(1 + x));
    
    % omega is the number of duplicates that can be made of the cell and
    % its hydrolases, given the local supply of nutrients
    omega = min([(kappa^2)*Ce/(CC*(1 + x)), ...
                 (kappa^2)*Ne/(Ni*(1 + x)), (kappa^2)*Pe/Pi]);
    
    % mu is our initial guess for the growth rate of mobile
    % cells when the relative density of exoenzymes is x 
    mu = (theta*phi*x/tau)*(1-exp(-omega))/...
            (1+phi*x-epsilon*exp(-omega));
        
    eta = mu;
    
    for dum = 1:sres
        
        mu = eta*(1 - exp(-mu*omega/eta));
       
        U = (theta*phi*x/tau) + epsilon*mu/(exp(mu*omega/eta)-1);
        
        eta = min(U, L)/(1 + phi*x);
    end
    
    if (isreal(mu) == 1) && (mu > 0)
        
        M_vector(i) = mu;
    end
end

[mu, eta, x] = find_best_autolytic...
    (Ci, Ni, Pi, Ce, Ne, Pe, kappa, tau, delta, L, epsilon, sres, xres);

figure(1)
plot(x_vector, M_vector)
hold on
plot(x,mu,'o')
xlabel('Relative density of exoenzymes')
ylabel('Apparent growth rate')