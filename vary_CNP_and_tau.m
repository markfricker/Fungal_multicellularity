%% set up parameters. These can be manually edited.
% res1 is the number of different values for the C:N ratio and recalcitrance
% that are tried
res1 = 200;

% xres determines the number of different values for x that are tried in finding
% the optimal solution. Final value found is accurate to xres^2.
xres = 200;

% sres is the number of iterations used in finding best solution for a given
% environment and a given rate of sythesis.
sres = 60;

% N_to_P set the N:P ratio of the resource.
N_to_P = 20;

% kappa is the radius of digestion, normalised to the radius of the cells/hypha,
% and sets the total amount of resource available.
kappa = 6;

% C_to_N_min and C_to_N_max define the limits of the C:N ratio in the resource.
% The number of intermediate steps is controlled by res1. The C:N ratio is
% plotted on the y-axis of the resultant map.
C_to_N_min = 5;
C_to_N_max = 300;

% tau_min and tau_max define the limits of the resource recalcitrance, where tau
% is the time in hours for hydrolases to digest their own mass. The number of
% intermediate steps is controlled by res1. The recalcitrance is plotted on the
% x-axis of the resultant map.
tau_min = 0.1;
tau_max = 100;

% Ci, Ni and Pi are the mass of carbon, nitrogen and phosphorous, respectively,
% needed per unit volume of growth, in g per ml. 
Ci = 0.165;
Ni = 0.032;
Pi = 0.005;

% the carbon use efficiency is the fraction of C for growth versus total C
% for respiration and growth
CUE = 0.5;

% Ct is the total C required, including respiratory C set by the carbon use
% efficiency
Ct = Ci/CUE;

% dry weight in grams per ml of the substrate.
density = 0.5;

% epsilon is the efficiency of recycling for senscent autolytic cells.
epsilon = 0.5;

% alpha is the additional mass of machinery needed for cell mobility, relative
% to the mass of essential metabolic machinery.
alpha = 0.02;

% beta is the mass of material in vesicles, relative to the mass of the rest of
% the fungus, including hydrolases.
beta = 0.1;

% lambda is the maximum rate of resource use per unit volume, in g per ml per
% hour.
lambda = 0.3;

% delta is the ratio of C that has to be digested to release each N, to reflect
% that N is embedded within C-rich polymers.
delta = 0;

%% set up arrays and derived parameters
% L is the maximum rate of resource use per unit volume, relative to Ct + Ni +
% Pi.
L = lambda/(Ct + Ni + Pi);

C_to_N_vector = zeros(res1, 1);
Recalcitrance = zeros(res1, 1);

Mu_immobile = zeros(res1);
Mu_motile = zeros(res1);
Mu_autolytic = zeros(res1);
Mu_fungal = zeros(res1);
Mu_cell = zeros(res1);

x_immobile = zeros(res1);
x_motile = zeros(res1);
x_autolytic = zeros(res1);
x_fungal = zeros(res1);
xC_fungal = zeros(res1);
xN_fungal = zeros(res1);
xP_fungal = zeros(res1);
x_cell = zeros(res1);

%% main program loop
for i = 1:res1
    
    C_to_N = exp(log(C_to_N_min) + ...
        (log(C_to_N_max) - log(C_to_N_min))*(i-1)/(res1-1));
    
    C_to_N_vector(i) = C_to_N;
    
    M_tot = 12*C_to_N*N_to_P + 14*N_to_P + 31;
    
    Ce = density*12*C_to_N*N_to_P/M_tot;
    Ne = density*14*N_to_P/M_tot;
    Pe = density*31/M_tot;
        
    for j = 1:res1
        
        tau = exp(log(tau_min) + ...
            (log(tau_max) - log(tau_min))*(j-1)/(res1-1));
        
        Recalcitrance(j) = tau;
        
        [Mu_cell(i,j), x_cell(i,j)] = find_best_cell(Ct, Ni, Pi, Ce, Ne, tau, delta, L);
        
        [Mu_immobile(i,j), Mu_immobile(i,j)] = find_best_immobile...
            (Ct, Ni, Pi, Ce, Ne, Pe, kappa, tau, delta, L, sres, xres);
        
        
        [Mu_motile(i,j), x_motile(i,j)] = find_best_motile...
            (Ct, Ni, Pi, Ce, Ne, Pe, kappa, tau, delta, ...
            L, alpha, sres, xres);
        
        [Mu_autolytic(i,j), ~, x_autolytic(i,j)] = find_best_autolytic...
            (Ct, Ni, Pi, Ce, Ne, Pe, kappa, tau, delta, ...
            L, epsilon, sres, xres);
        
        [Mu_fungal(i,j), x_fungal(i,j), x_fungal(i,j), ...
            xN_fungal(i,j), xP_fungal(i,j)] = find_best_fungi(Ct, Ni, ...
            Pi, Ce, Ne, Pe, kappa, tau, delta, L, beta, sres, 3*xres);
        
    end
    
    percent_finished = 100*i/res1
end

clear tau
clear M_tot
clear percent_finished

