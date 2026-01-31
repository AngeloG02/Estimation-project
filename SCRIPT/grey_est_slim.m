% function [theta_est,std_theta] = grey_est_slim (sigma_q,sigma_ax)

clc; clear;close all


main_quad_ANTX % inizializza i parametri che servono al modello


%%  INIZIALIZATION
f_min = 0.85;
f_max = 1.9;

% CAMBIARE IN BASE ALLA CONFIGURAZIONE MIGLIORE
sigma_q = 1; % [deg/s]
sigma_ax = 0.5; % [m/s^2]


% Noise
noise.Enabler = 1;
noise.pos_stand_dev = noise.Enabler * 0.0011;
noise.vel_stand_dev = noise.Enabler * 0.01;
noise.attitude_stand_dev = noise.Enabler * deg2rad(0.33);
noise.acc_stand_dev = noise.Enabler * sigma_ax;
noise.ang_rate_stand_dev = noise.Enabler * deg2rad(sigma_q);

assignin('base', 'noise',noise);

T = 0;
aux = {};
theta_0 = [-0.5; -0.1; 0.5; -20; -1; 1];
init_sys = idgrey(@System_matrix, theta_0, 'c', aux, T);
data.InputName = 'Pitch Moment';
data.InputUnit = '-';
data.OutputName = {'Pitch rate', 'Longitudinal acceleration'};
data.OutputUnit = {'rad/s', 'm/s^2'};
data.TimeUnit = 's';

opt = greyestOptions;
opt.InitialState = 'backcast';
opt.Focus = 'simulation';
opt.Display = 'on';  
opt.EnforceStability = false;
opt.SearchMethod = 'auto';
opt.SearchOptions.MaxIterations = 100;

theta_true = [0.1068 0.1192 -5.9755 -2.6478 -10.1647 450.71];
for i = 1:length(theta_true)
    if theta_true(i) < 0
        init_sys.Structure.Parameters(1).Minimum(i) = theta_true(i)* 2;
        init_sys.Structure.Parameters(1).Maximum(i) = theta_true(i) /2;
    else
        init_sys.Structure.Parameters(1).Minimum(i) = theta_true(i) /2;
        init_sys.Structure.Parameters(1).Maximum(i) = theta_true(i) *2;
    end
end


function [A, B, C ,D]= System_matrix(theta,T)
    theta = theta(:);
    X_u = theta(1);
    X_q = theta(2);
    M_u = theta(3);
    M_q = theta(4);
    X_d = theta(5);
    M_d = theta(6);
    g = 9.81;
    
    A = [X_u, X_q, -g;
         M_u, M_q,  0;
           0,   1,  0];
    B = [X_d; M_d; 0];
    C = [0, 1, 0; X_u, X_q, 0];
    D = [0 ; X_d];
end




%% STEP 1: DATA ACQUISITION (Multiple Realizations)
N_realizations = 15;
theta_est = zeros(6, N_realizations);
std_theta = zeros(6, N_realizations);

% ===== CREA ARRAY DI SimulationInput =====
simIn(1:N_realizations) = Simulink.SimulationInput('Simulator_Single_Axis');

for i = 1:N_realizations
    % Imposta i seed per ogni simulazione
    simIn(i) = simIn(i).setVariable('Seed_pos', i);
    simIn(i) = simIn(i).setVariable('Seed_vel', i + 1);
    simIn(i) = simIn(i).setVariable('Seed_theta', i + 2);
    simIn(i) = simIn(i).setVariable('Seed_q', i + 3);
    simIn(i) = simIn(i).setVariable('Seed_ax', i + 4);
end

% ===== ESEGUI SIMULAZIONI IN PARALLELO =====
simOut = parsim(simIn, 'ShowProgress', 'on');

% ===== PROCESSA I RISULTATI =====
parfor realization = 1:N_realizations
    fprintf('  Realization %2d/%d: Processing... ', realization, N_realizations);
    
    % Estrai i risultati dalla simulazione
    out = simOut(realization);
    
    t = out.q.Time;
    u = out.Mtot.Data;
    q = out.q.Data;
    ax = out.ax.Data;
    
    Ts = t(2) - t(1);
    
    y = [q, ax];
    N = length(y);
    Y = fft(y);
    U = fft(u);
    
    freq_Hz = (0:N-1)' / (N * Ts);
    idx = (freq_Hz >= f_min) & (freq_Hz <= f_max);
    
    Y_filtered = Y(idx, :);
    U_filtered = U(idx);
    W_filtered = 2*pi * freq_Hz(idx);
    
    data = iddata(Y_filtered, U_filtered, 0, 'Frequency', W_filtered);
    [sys, ~] = greyest(data, init_sys, opt);
    
    theta_est(:, realization) = getpvec(sys);
    cov_theta = getcov(sys);
    std_theta(:, realization) = sqrt(diag(cov_theta));
    
    fprintf('Done\n');
end

% end

