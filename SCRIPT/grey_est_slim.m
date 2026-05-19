% function [theta_est,std_theta] = grey_est_slim (sigma_q,sigma_ax)

% clc; clear;close all


main_quad_ANTX % inizializza i parametri che servono al modello


%%  INIZIALIZATION
f_min = f_min_sx;
f_max = f_max_sx;

% CAMBIARE IN BASE ALLA CONFIGURAZIONE MIGLIORE
sigma_q = sigma_q_sx;     % [deg/s] % Caso SINISTRA
sigma_ax = sigma_ax_sx;   % [m/s^2] % Caso SINISTRA


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
theta_0 = [-0.001; 0.001; -5.8; -4.7; -10; 120];
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
N_realizations = 50;
theta_est = zeros(6, N_realizations);
std_theta = zeros(6, N_realizations);
theta_mc = zeros(6, N_realizations); % 6 parametri
poles_mc = zeros(3, N_realizations); % 3 poli per la dinamica longitudinale
zeros_mc = cell(1, N_realizations);  % cell array perché il numero di zeri può variare
sys_mc   = cell(1, N_realizations);  % salviamo tutti i modelli per il Bode

% ===== CREA ARRAY DI SimulationInput =====
simIn(1:N_realizations) = Simulink.SimulationInput('Simulator_Single_Axis');

for i = 1:N_realizations
    % Imposta i seed per ogni simulazione
    simIn(i) = simIn(i).setVariable('Seed_pos', i);
    simIn(i) = simIn(i).setVariable('Seed_vel', i + N_realizations);
    simIn(i) = simIn(i).setVariable('Seed_theta', i + 2*N_realizations);
    simIn(i) = simIn(i).setVariable('Seed_q', i + 3*N_realizations);
    simIn(i) = simIn(i).setVariable('Seed_ax', i + 4*N_realizations);
end

% ===== ESEGUI SIMULAZIONI IN PARALLELO =====
simOut = parsim(simIn, 'ShowProgress', 'on', 'TransferBaseWorkspaceVariables', 'on');

% ===== PROCESSA I RISULTATI =====
parfor realization = 1:N_realizations
    fprintf('  Realization %2d/%d: Processing... ', realization, N_realizations);
    
    % Estrai i risultati dalla simulazione
    out = simOut(realization);
    
    t = out.q.Time;
    u = out.Mtot.Data;
    q = out.q.Data;
    ax = out.ax.Data;
    
    % Rimozione della media
    u = u - mean(u);
    q = q - mean(q);
    ax = ax - mean(ax);

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
    
    data = iddata(Y_filtered, U_filtered, 0.004, 'Frequency', W_filtered);
    [sys, ~] = greyest(data, init_sys, opt);
    
    % salvataggio incertezza
    theta_est(:, realization) = getpvec(sys);
    cov_theta = getcov(sys);
    std_theta(:, realization) = sqrt(diag(cov_theta));
    
    % salvataggio per i grafici montecarlo
    theta_mc(:, realization) = getpvec(sys);
    poles_mc(:, realization) = pole(sys);
    zeros_mc{realization}    = tzero(sys);
    sys_mc{realization}      = sys;

    fprintf('Done\n');
end

% end
%% ========================================================================
%  TASK 3.2 - STATISTICAL ANALYSIS (MONTE CARLO)
% =========================================================================

% 1. ISTOGRAMMI DEI PARAMETRI STIMATI
figure('Name', 'Monte Carlo: Parameter Histograms');
param_names = {'X_u', 'X_q', 'M_u', 'M_q', 'X_\delta', 'M_\delta'};
for j = 1:6
    subplot(2, 3, j);
    % 'pdf' per normalizzare l'istogramma come una densità di probabilità
    histogram(theta_mc(j, :), 15, 'Normalization', 'pdf', 'FaceColor', [0.2 0.6 0.8]);
    title(sprintf('Distribuzione di %s', param_names{j}));
    xlabel('Valore'); ylabel('Densità di probabilità');
    grid on;
end

% 2. MAPPA DEI POLI E DEGLI ZERI (Dispersione nel piano complesso)
figure('Name', 'Monte Carlo: Poles and Zeros Dispersion');
hold on; grid on;
% Plotta tutti i poli in blu
plot(real(poles_mc(:)), imag(poles_mc(:)), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5);
% Plotta tutti gli zeri in rosso
for j = 1:N_realizations
    z = zeros_mc{j};
    if ~isempty(z)
        plot(real(z), imag(z), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
    end
end
xline(0, 'k-', 'LineWidth', 1.5); % Asse immaginario (limite di stabilità)
xlabel('Asse Reale'); ylabel('Asse Immaginario');
title('Dispersione di Poli (x) e Zeri (o)');
legend('Poli', 'Zeri', 'Location', 'best');

% 3. DISPERSIONE DELLA RISPOSTA IN FREQUENZA (FRF)
figure('Name', 'Monte Carlo: FRF Dispersion');
% L'operatore {:} espande la cell array passando tutti i 50 modelli a bode
bode(sys_mc{:});
grid on;
title('Dispersione della FRF (50 Realizzazioni Monte Carlo)');
