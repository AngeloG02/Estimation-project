
clc; clear;close all

%% cento
main_quad_ANTX % inizializza i parametri, non runna il modello 
sigma_q_cent = 1; % [deg/s]
sigma_ax_cent = 0.5; % [m/s^2]
[f_min_cent,f_max_cent] = best_freq_calculator(sigma_q_cent,sigma_ax_cent);

main_quad_ANTX % inizializza i parametri, non runna il modello 
[theta_est_cent,std_theta_cent,cov_cent] = grey_est(sigma_q_cent,sigma_ax_cent,f_min_cent,f_max_cent);

% main_quad_ANTX
% sigma_q_cent = 1; % [deg/s]
% sigma_ax_cent = 0.5; % [m/s^2]
% [theta_est_cent_manuale,std_theta_cent_manuale] = grey_est_manual(sigma_q_cent,sigma_ax_cent);

%% sx
main_quad_ANTX % inizializza i parametri, non runna il modello 
sigma_q_sx= 1.3; % [deg/s]
sigma_ax_sx = 0.35; % [m/s^2]
[f_min_sx,f_max_sx] = best_freq_calculator(sigma_q_sx,sigma_ax_sx);

main_quad_ANTX % inizializza i parametri, non runna il modello (ho commentato)
sigma_q_sx= 1.3; % [deg/s]
sigma_ax_sx = 0.35; % [m/s^2]
[theta_est_sx,std_theta_sx,cov_sx] = grey_est(sigma_q_sx,sigma_ax_sx,f_min_sx,f_max_sx);

%% dx
main_quad_ANTX % inizializza i parametri, non runna il modello 
sigma_q_dx  = 0.7; % [deg/s]
sigma_ax_dx = 0.65; % [m/s^2]
[f_min_dx,f_max_dx] = best_freq_calculator(sigma_q_dx,sigma_ax_dx);

main_quad_ANTX % inizializza i parametri, non runna il modello (ho commentato)
sigma_q_dx  = 0.7; % [deg/s]
sigma_ax_dx = 0.65; % [m/s^2]
[theta_est_dx,std_theta_dx,cov_dx] = grey_est(sigma_q_dx,sigma_ax_dx,f_min_dx,f_max_dx);




% Scegli la configurazione migliore e usala nel grey_est_slim


%% ========================================================================
%  CONFRONTO INCERTEZZE (TASK 2.2) - I 3 CRITERI DI EXPERIMENT DESIGN
% =========================================================================
casi = {'CENTRALE', 'SINISTRA', 'DESTRA'};

% --- 1. CRITERIO DELLA TRACCIA (A-optimality): min Tr[cov] ---
% Misura l'incertezza media globale (Somma delle varianze)
Tr_cent = trace(cov_cent);
Tr_sx   = trace(cov_sx);
Tr_dx   = trace(cov_dx);

% --- 2. CRITERIO DEL DETERMINANTE (D-optimality): min Det[cov] ---
% Misura il volume dell'ellissoide di incertezza nello spazio a 6 dimensioni
Det_cent = det(cov_cent);
Det_sx   = det(cov_sx);
Det_dx   = det(cov_dx);

% --- 3. CRITERIO DELL'AUTOVALORE MASSIMO (E-optimality): min max(eig[cov]) ---
% Misura l'incertezza nella direzione "peggiore" possibile
Eig_cent = max(eig(cov_cent));
Eig_sx   = max(eig(cov_sx));
Eig_dx   = max(eig(cov_dx));

fprintf('\n\n================ CONFRONTO CRITERI DI EXPERIMENT DESIGN ================\n');

fprintf('\n1. CRITERIO DELLA TRACCIA (minima incertezza media globale)\n');
fprintf('   Centro: %10.4e  |  Sinistra: %10.4e  |  Destra: %10.4e\n', Tr_cent, Tr_sx, Tr_dx);
[~, b1] = min([Tr_cent, Tr_sx, Tr_dx]);
fprintf('   >>> VINCITORE TRACCIA: Caso %s\n', casi{b1});

fprintf('\n2. CRITERIO DEL DETERMINANTE (minimo volume dell''ellissoide di incertezza)\n');
fprintf('   Centro: %10.4e  |  Sinistra: %10.4e  |  Destra: %10.4e\n', Det_cent, Det_sx, Det_dx);
[~, b2] = min([Det_cent, Det_sx, Det_dx]);
fprintf('   >>> VINCITORE DETERMINANTE: Caso %s\n', casi{b2});

fprintf('\n3. CRITERIO DELL''AUTOVALORE MASSIMO (minima incertezza nel caso peggiore)\n');
fprintf('   Centro: %10.4e  |  Sinistra: %10.4e  |  Destra: %10.4e\n', Eig_cent, Eig_sx, Eig_dx);
[~, b3] = min([Eig_cent, Eig_sx, Eig_dx]);
fprintf('   >>> VINCITORE AUTOVALORE MAX: Caso %s\n', casi{b3});
fprintf('\n========================================================================\n');




%%
%%%%%%%%% MONTECARLO %%%%%%%%%%%%%%%%%%%%%%
grey_est_slim 

% fare statistica su theta_est
