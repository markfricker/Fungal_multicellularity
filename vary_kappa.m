% res1 is the number of different values for kappa that we try
res1 = 200;

% xres is the number of iterations used in finding best solution for 
% a given environment and a given rate of sythesis
xres = 200;

% sres is the number of iterations used in finding best solution for 
% a given environment and a given rate of sythesis
sres = 60;

% k_min and k_max are the minimum and maximum values for the resource
% acquisition length that we try
k_min = 1;
k_max = 8;

% Ci is the mass of carbon needed per unit volume of growth, in g per ml
Ci = 0.33;
Ni = 0.032;
Pi = 0.005;

C_to_P = 2000;
N_to_P = 10;

% tau is the time in hours for hydrolases to digest their own mass
tau = 40;

% dry weight in grams per ml of the substrate
density = 0.5;

M_tot = 12*C_to_P + 14*N_to_P + 31;
Ce = density*12*C_to_P/M_tot;
Ne = density*14*N_to_P/M_tot;
Pe = density*31/M_tot;

% delta is the minimal fraction of C that must be digested in order to 
% digest the available N
delta = 0;

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

% phi is the correction term
phi = (Ci + Ni)/(Ci + Ni + Pi);

kappa_vector = zeros(res1, 1);

% apparent growth rates for each category of organism
Mu_immobile = zeros(res1, 1);
Mu_motile = zeros(res1, 1);
Mu_autolytic = zeros(res1, 1);
Mu_fungal = zeros(res1, 1);

% rate of digestion per unit volume for each category of organism
D_immobile = zeros(res1, 1);
D_motile = zeros(res1, 1);
D_autolytic = zeros(res1, 1);
D_fungal = zeros(res1, 1);
DC_fungal = zeros(res1, 1);
DP_fungal = zeros(res1, 1);
DN_fungal = zeros(res1, 1);

for i = 1:res1
    
    kappa = k_min + (k_max - k_min)*(i-1)/(res1-1);
    kappa_vector(i) = kappa;
    
    [M, x] = find_best_immobile...
    (Ci, Ni, Pi, Ce, Ne, Pe, kappa, tau, delta, L, sres, xres);
     
    Mu_immobile(i) = M;
    D_immobile(i) = phi*x/tau;
    
    [M, x] = find_best_motile(Ci, Ni, Pi, Ce, Ne, Pe, ...
        kappa, tau, delta, L, alpha, sres, xres);
    
    Mu_motile(i) = M;
    D_motile(i) = phi*x/tau;
    
    [M, ~, x] = find_best_autolytic(Ci, Ni, Pi, Ce, Ne, Pe, ...
        kappa, tau, delta, L, epsilon, sres, xres);

    Mu_autolytic(i) = M;
    D_autolytic(i) = phi*x/tau;
    
    [M, x, xC, xN, xP] = find_best_fungi...
        (Ci, Ni, Pi, Ce, Ne, Pe, kappa, tau, delta,L,beta,sres,3*xres);
    
    Mu_fungal(i) = M;
    DC_fungal(i) = phi*xC/tau;
    DN_fungal(i) = phi*xN/tau;
    DP_fungal(i) = phi*xP/tau;
    D_fungal(i) = phi*x/tau;
    
    clc
    percent_finished = 100*i/res1
end

figure
plot(kappa_vector, Mu_motile*24, 'c', 'LineWidth', 2)
hold on
plot(kappa_vector, Mu_autolytic*24, 'g', 'LineWidth', 2)
plot(kappa_vector, Mu_fungal*24, 'm', 'LineWidth', 2)
plot(kappa_vector, Mu_immobile*24, 'k--', 'LineWidth', 2)
xlabel('Digestion Range \kappa')
ylabel('Apparent growth rate, day^{-1}')
legend('Motile', 'Autolytic', 'Fungi', 'Immobile')


figure
plot(kappa_vector, D_motile*(Ci+Ni+Pi)*24, 'c', 'LineWidth', 2)
hold on
plot(kappa_vector, D_autolytic*(Ci+Ni+Pi)*24, 'g', 'LineWidth', 2)
plot(kappa_vector, D_fungal*(Ci+Ni+Pi)*24, 'm', 'LineWidth', 2)
plot(kappa_vector, D_immobile*(Ci+Ni+Pi)*24, 'k--', 'LineWidth', 2)
xlabel('Digestion Range \kappa')
ylabel('Digestion Rate, g ml^{-1} day^{-1}')
legend('Motile', 'Autolytic', 'Fungi', 'Immobile')

