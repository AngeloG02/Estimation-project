
clc; clear;close all
main_quad_ANTX % inizializza i parametri, non runna il modello (ho commentato)

sigma_q_cent = 1; % [deg/s]
sigma_ax_cent = 0.5; % [m/s^2]
[theta_est_cent,std_theta_cent] = grey_est(sigma_q_cent,sigma_ax_cent);

main_quad_ANTX % inizializza i parametri, non runna il modello (ho commentato)
sigma_q_sx= 1.3; % [deg/s]
sigma_ax_sx = 0.35; % [m/s^2]
[theta_est_sx,std_theta_sx] = grey_est(sigma_q_sx,sigma_ax_sx);

main_quad_ANTX % inizializza i parametri, non runna il modello (ho commentato)
sigma_q_dx  = 0.7; % [deg/s]
sigma_ax_dx = 0.65; % [m/s^2]
[theta_est_dx,std_theta_dx] = grey_est(sigma_q_dx,sigma_ax_dx);




% Scegli la configurazione migliore e usala nel grey_est_slim




%%
%%%%%%%%% MONTECARLO %%%%%%%%%%%%%%%%%%%%%%
grey_est_slim 

% fare statistica su theta_est
