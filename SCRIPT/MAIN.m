clc
clear all
close all



%% STEP 1: DATA ACQUISITION (Multiple Realizations)
fprintf('\nStep 1: Acquiring data from Simulink runs...\n');

N_realizations = 1;      % non mettere più di una (per ora) 

% Preallocate storage for concatenated signals
u_combined = [];
q_combined = [];
ax_combined = [];
time_offset = 0;
t_combined = [];
fs = [];  % Will be extracted from first run

% Loop over realizations with different noise seeds
for realization = 1:N_realizations
    fprintf('  Realization %2d/%d: ', realization, N_realizations);
    
    % Run simulation
    main_quad_ANTX;
    % Extract signals from this Simulink run
    t_single = out.q.Time;
    u_single = out.Mtot.Data;      % Input: Momento totale
    q_single = out.q.Data;         % Output 1: Pitch rate
    ax_single = out.ax.Data;       % Output 2: Accel. longitudinale
    
    % Get sampling info from first run
    if realization == 1
        Ts = t_single(2) - t_single(1);
        fs = 1/Ts;
        samples_per_run = length(t_single);
        fprintf('(fs=%.1f Hz, %d samples per run)', fs, samples_per_run);
    end
    
    % Adjust time vector for concatenation (reset to 0 each run)
    t_single_adjusted = t_single + time_offset;
    
    % Concatenate signals from this realization
    u_combined = [u_combined; u_single];
    q_combined = [q_combined; q_single];
    ax_combined = [ax_combined; ax_single];
    t_combined = [t_combined; t_single_adjusted];
    % CELLA per ogni esperimento (NON concatenare)
    u_experiments{realization} = u_single;
    q_experiments{realization} = q_single;
    ax_experiments{realization} = ax_single;
    t_experiments{realization} = t_single_adjusted;

    
    % Update time offset for next realization
    time_offset = time_offset + t_single(end) + Ts;
    fprintf('\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MANUAL APPROACH %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL SYSTEM FRF

G = @(theta,freq_vect)Model_system_FRF(theta,freq_vect);
%% REAL SYSTEM FRF

[H_q_trusted,H_ax_trusted,freq_vect] = FRF_real_system_estimation(u_combined,q_combined,ax_combined,t_combined,fs);

G_m = [H_q_trusted,H_ax_trusted];
 
%% OTTIMIZZATION
% theta_0 = [-0.5; -0.1; 0.5; -20; -1; 1];

% OPTIMIZATION_WITH_FMINCON
% [theta_hat_oe, J_hat_oe] = outputError_FD(G_m, freq_vect, G, theta_0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AUTOMATIC APPROACH %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% =============================================
%  STEP 0: Generiamo i dati
%  =============================================
% Suppongo che tu abbia già: u_experiments, q_experiments, ax_experiments, N_realizations, Ts

% Costruisci i dati multi-esperimento
aux = {};
T = 0;

% for i = 1:N_realizations
%     y_cell{i} = [q_experiments{i}, ax_experiments{i}];
%     u_cell{i} = u_experiments{i};
% end
% y = cell2mat(y_cell);
% u = cell2mat(u_cell);

y = [q_combined, ax_combined];
u = u_combined;

N = length(y);
Y = fft(y);
U = fft(u);

% Procedure: Isolate dynamic response by removing the mean
    % Rationale: Spectral analysis focuses on fluctuations. Removing the DC
    % offset prevents a large, distorting spike at 0 Hz in the spectrum.
    u = u_combined - mean(u_combined);
    q = q_combined - mean(q_combined);
    ax = ax_combined - mean(ax_combined);

%-->Step#3:Spectral Estimation (Welch's Method)

    % Apply the Welch method to compute smooth estimates of the auto- and
    % cross-power spectral densities (PSDs). This robustly handles the high
    % variance of a simple periodogram.
    fprintf('Step 3: Performing spectral estimation using Welch''s method...\n');
    % Define Welch method parameters
    M = 2048;               % Segment length (a power of 2 for FFT efficiency)
    noverlap = M / 2;       % 50% overlap between segments
    win = hanning(M);       % Hanning window to reduce spectral leakage

    % Compute the required spectra using MATLAB's built-in functions
    % 'pwelch' and 'cpsd' implement the entire 4-step Welch procedure:
    % (Segment -> Window -> Transform -> Average)

    % Input Autospectrum (Suu)
    [Suu, f] = pwelch(u, win, noverlap, [], fs);
    % Cross-spectrum between input (u) and pitch rate (q)
    [Suq, ~] = cpsd(u, q, win, noverlap, [], fs);
    % Cross-spectrum between input (u) and longitudinal acceleration (ax)
    [Suax, ~] = cpsd(u, ax, win, noverlap, [], fs);

%  %-->Step#4:FRF Computation 
%  % Compute the non-parametric FRF estimate for each output channel using
%  % the fundamental spectral relationship: G(f) = S_uy(f) / S_uu(f).
% 
% 
% 
%-->Step#5:Quality Assesment   
    % Compute the coherence function for each input-output pair to
    % quantitatively assess the reliability of the FRF estimate.
    fprintf('Step 5: Assessing FRF quality with the coherence function...\n');
    % To compute coherence, we first need the output autospectra (Sqq, Saxax)
    [Sqq, ~] = pwelch(q, win, noverlap, [], fs);
    [Saxax, ~] = pwelch(ax, win, noverlap, [], fs);
    % Coherence formula: gamma^2 = |S_uy|^2 / (S_uu * S_yy)
    gamma_sq_q = (abs(Suq).^2) ./ (Suu .* Sqq);
    gamma_sq_ax = (abs(Suax).^2) ./ (Suu .* Saxax);


% Aggiungi linee verticali per la banda di interesse (MODIFICARE MANUALMENTE)
f_min = 0.85;   % ← MODIFICA PER PITCH RATE
f_max = 1.9;   % ← MODIFICA PER PITCH RATE


%-->Step#6:Vizualization   
    % Plot the results for analysis.
    fprintf('Step 6: Generating plots...\n');

    % Plot 1: Estimated FRFs (Bode Plot)
    figure
    % Pitch Rate Channel
    subplot(2,2,1);
    semilogx(f, 20*log10(abs(Suq)));
    grid on; title('CSPD q '); ylabel('Magnitude [dB]');
    xline(f_min, 'r--', 'LineWidth', 2.5, 'Label', sprintf('f_{min}=%.2f Hz', f_min));
    xline(f_max, 'g--', 'LineWidth', 2.5, 'Label', sprintf('f_{max}=%.2f Hz', f_max));

    % Longitudinal Acceleration Channel
    subplot(2,2,2);
    semilogx(f, 20*log10(abs(Suax)));
    grid on; title('CSPD q ');
    xline(f_min, 'r--', 'LineWidth', 2.5, 'Label', sprintf('f_{min}=%.2f Hz', f_min));
    xline(f_max, 'g--', 'LineWidth', 2.5, 'Label', sprintf('f_{max}=%.2f Hz', f_max));

    % Pitch Rate coherence
    subplot(2,2,3);
    semilogx(f, gamma_sq_q);
    hold on;
    plot([f(1), f(end)], [0.7, 0.7], 'r--', 'LineWidth', 1.5); % Threshold line
    grid on; title('Coherence: Input (u) to Pitch Rate (q)');
    ylim([0 1.05]); ylabel('Coherence \gamma^2'); legend('Coherence', 'Trust Threshold');
    xline(f_min, 'r--', 'LineWidth', 2.5, 'Label', sprintf('f_{min}=%.2f Hz', f_min));
    xline(f_max, 'g--', 'LineWidth', 2.5, 'Label', sprintf('f_{max}=%.2f Hz', f_max));

    % Longitudinal Acceleration coherence
    subplot(2,2,4);
    semilogx(f, gamma_sq_ax);
    hold on;
    plot([f(1), f(end)], [0.7, 0.7], 'r--', 'LineWidth', 1.5); % Threshold line
    grid on; title('Coherence: Input (u) to Longitudinal Accel. (ax)');
    ylim([0 1.05]); xlabel('Frequency [Hz]'); ylabel('Coherence \gamma^2');
    xline(f_min, 'r--', 'LineWidth', 2.5, 'Label', sprintf('f_{min}=%.2f Hz', f_min));
    xline(f_max, 'g--', 'LineWidth', 2.5, 'Label', sprintf('f_{max}=%.2f Hz', f_max));
 %%

idx_f_trusted = (f >= f_min) & (f <= f_max);
freq_trusted = f(idx_f_trusted);


%% =============================================
%  STEP 5: CREA MODELLO INIZIALE
%  =============================================
theta_0 = [-0.5; -0.1; 0.5; -20; -1; 1];
init_sys = idgrey(@System_matrix, theta_0, 'c', aux, T);

% Controlla stabilità
eig_init = eig(init_sys.A);
fprintf('\n========== INITIAL MODEL CHECK ==========\n');
fprintf('Eigenvalues of A:\n');
disp(eig_init);

if all(real(eig_init) < 0)
    fprintf('✓ Model is STABLE\n');
else
    fprintf('✗ WARNING: Model is UNSTABLE!\n');
    fprintf('  Adjust theta_0 values.\n');
end

% Vettore di frequenze (in Hz, non rad/s per semplicità)
freq_Hz = (0:N-1)' / (N * Ts);

% =============================================
%  STEP 2: SELEZIONA INTERVALLO DI FREQUENZE
%  =============================================
% Definisci il range di interesse (in Hz)
f_min = freq_trusted(1);   % frequenza minima (Hz)
f_max = freq_trusted(end);     % frequenza massima (Hz)

% Trova gli indici dentro l'intervallo
idx = (freq_Hz >= f_min) & (freq_Hz <= f_max);

% Estrai solo i dati nell'intervallo
Y_filtered = Y(idx, :);  % Mantiene entrambe le colonne
U_filtered = U(idx);
W_filtered = 2*pi * freq_Hz(idx);  % Frequenze in rad/s

fprintf('Frequency range: %.2f - %.2f Hz\n', f_min, f_max);
fprintf('Number of frequency points: %d\n', sum(idx));

% =============================================
%  STEP 3: Crea iddata
%  =============================================
data = iddata(Y_filtered, U_filtered,0, 'Frequency', W_filtered);

data.InputName = 'Pitch Moment';
data.InputUnit = '-';
data.OutputName = {'Pitch rate', 'Longitudinal acceleration'};
data.OutputUnit = {'rad/s', 'm/s^2'};
data.TimeUnit = 's';




%% =============================================
%  STEP 6: CONFIGURA OPZIONI
%  =============================================
opt = greyestOptions;
opt.InitialState = 'backcast';
opt.Focus = 'simulation';
opt.Display = 'on';
opt.EnforceStability = false;
opt.SearchMethod = 'auto';
opt.SearchOptions.MaxIterations = 200;

% assegna un range in cui cercare i parametri
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

%% =============================================
%  STEP 7: ESEGUI GREYEST CON DATI FREQUENZIALI
%  =============================================
fprintf('\n========== RUNNING GREYEST (FREQUENCY DOMAIN) ==========\n');
[sys, x0] = greyest(data, init_sys, opt);

fprintf('\n========== ESTIMATION COMPLETED ==========\n');

%% =============================================
%  STEP 8: ANALIZZA RISULTATI
%  =============================================
fprintf('\n========== ESTIMATED PARAMETERS ==========\n');
theta_est = getpvec(sys);

param_names = {'X_u', 'X_q', 'M_u', 'M_q', 'X_d', 'M_d'};
for i = 1:length(theta_est)
    fprintf('%s = %.6f\n', param_names{i}, theta_est(i));
end

for i = 1:length(theta_true)
    error_pct = abs(theta_est(i) - theta_true(i)) / abs(theta_true(i)) * 100;
    fprintf('theta(%d)           % -15.6f % -15.6f % -15.2f\n', ...
        i, theta_true(i), theta_est(i), error_pct);
end

eig_est = eig(sys.A);
disp(eig_est);


fprintf('\n========== FIT REPORT ==========\n');
disp(sys.Report.Fit);

fprintf('\n========== COVARIANCE & UNCERTAINTY ==========\n');
try
    [cov_theta] = getcov(sys);
    std_theta = sqrt(diag(cov_theta));
    fprintf('Parameter uncertainties (std dev):\n');
    for i = 1:length(theta_est)
        fprintf('%s = %.6f ± %.6f\n', param_names{i}, theta_est(i), std_theta(i));
    end
catch
    fprintf('Covariance matrix not available\n');
end

%% T1.2 - Step 2: Confidence Intervals

z_1sigma = 1.000;   % 68% CI
z_2sigma = 1.960;   % 95% CI (quello più importante)
z_3sigma = 2.576;   % 99% CI

fprintf('\n95%% Confidence Intervals:\n');
for i = 1:length(theta_est)
    lower = theta_est(i) - z_2sigma * std_theta(i);
    upper = theta_est(i) + z_2sigma * std_theta(i);
    fprintf('%s: [%.6f, %.6f]\n', param_names{i}, lower, upper);
end

%% T1.2 - Step 3: Parameter Correlation

% Correlazione: ρ_ij = cov(θ_i, θ_j) / (σ_i * σ_j)
rho = cov_theta ./ (std_theta * std_theta');

fprintf('\nParameter Correlation Matrix:\n');
disp(rho);

% Identifica parametri altamente correlati
for i = 1:length(theta_est)
    for j = i+1:length(theta_est)
        if abs(rho(i,j)) > 0.9
            fprintf('⚠️ HIGH CORRELATION: %s ↔ %s (ρ = %.3f)\n', ...
                param_names{i}, param_names{j}, rho(i,j));
        end
    end
end
%% T1.2 - Step 4: Identifiability Assessment

relative_std = 100 * std_theta ./ abs(theta_est);
[~, idx_ranked] = sort(relative_std, 'descend');

fprintf('\nParameter Identifiability Ranking:\n');
for rank = 1:length(theta_est)
    i = idx_ranked(rank);
    if relative_std(i) < 5
        status = 'Excellent';
    elseif relative_std(i) < 10
        status = 'Good';
    elseif relative_std(i) < 25
        status = 'Moderate';
    else
        status = 'Poor';
    end
    fprintf('%d. %s: %.2f%% → %s\n', rank, ...
        param_names{i}, relative_std(i), status);
end
%% T1.2 - Step 5: Validation Against True Values

z_2sigma = 1.96;
theta_true = [0.1068 0.1192 -5.9755 -2.6478 -10.1647 450.71];
coverage_count = 0;

fprintf('\n95%% CI Coverage of True Parameters:\n');
for i = 1:length(theta_est)
    lower = theta_est(i) - z_2sigma * std_theta(i);
    upper = theta_est(i) + z_2sigma * std_theta(i);
    
    if (theta_true(i) >= lower) && (theta_true(i) <= upper)
        fprintf('✓ %s: TRUE value INSIDE CI\n', param_names{i});
        coverage_count = coverage_count + 1;
    else
        fprintf('✗ %s: TRUE value OUTSIDE CI\n', param_names{i});
    end
end

coverage_pct = 100 * coverage_count / length(theta_est);
fprintf('\nCoverage: %d/%d (%.0f%%)\n', coverage_count, length(theta_est), coverage_pct);

if coverage_pct >= 90
    fprintf('✓ Uncertainty correctly estimated\n');
elseif coverage_pct >= 60
    fprintf('⚠️ Some underestimation of uncertainty\n');
else
    fprintf('✗ Significant underestimation\n');
end

%% =============================================
%  STEP 9: BODE PLOT DEL MODELLO IDENTIFICATO
%  =============================================
figure;
bode(sys);
grid on;
title('Bode Plot of Identified Grey-Box Model (Frequency Domain Estimation)', 'FontSize', 14);

%% =============================================
%  STEP 10: CONFRONTO MODELLO vs DATI FREQUENZIALI
%  =============================================
figure;
opt_comp = compareOptions('InitialCondition', 'zero');
compare(data, sys, opt_comp);


%% 

function [A, B, C ,D]= System_matrix(theta,T)
% restituisce le due FRF (uscita1 = q, uscita2 = a_x)

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

    B = [X_d;
         M_d;
           0];

    C = [0,   1, 0;      % prima uscita (q)
         X_u, X_q, 0];   % seconda uscita (a_x)

    D = [0 ; X_d];        % [2x1] feedthrough


   
end

%% T1.3 validation with 3211 sequence

% 3211 SEQUENCE GENERATION (τ = 2 s, A = 0.1)

Ts   = 0.01;          % sampling time [s]  (modifica se necessario)
Tend = 18;            % total duration [s]
t_3211 = (0:Ts:Tend)'; 

A = 0.1;              % amplitude
tau = 2;              % step duration [s]

u_3211 = zeros(size(t_3211));

% pattern 3-2-1-1:
%  0–2 s   : +A
%  2–4 s   :  0
%  4–8 s   : +A
%  8–12 s  : -A
% 12–14 s  : +A
% 14–16 s  : -A
% 16–18 s  :  0

u_3211(t_3211 >= 0  & t_3211 < tau) =  0;   % 0–4 s   (3 e 1)
u_3211(t_3211 >= tau & t_3211 < 4*tau) =  A;   % 4–8 s   (2)
u_3211(t_3211 >= 4*tau & t_3211 < 6*tau) = -A;   % 8–12 s  (1)
u_3211(t_3211 >= 6*tau & t_3211 < 7*tau) =  A;   % 12–14 s (1)
u_3211(t_3211 >= 7*tau & t_3211 < 8*tau) = -A;   % 14–16 s
u_3211(t_3211 >= 8*tau) = 0;                    % 16–18 s

% (facoltativo) visualizza
figure; stairs(t_3211, u_3211, 'LineWidth', 2); grid on;
xlabel('Time [s]'); ylabel('Input amplitude'); title('3211 Sequence');

clear ExcitationM
ExcitationM(:,1)= t_3211; 
ExcitationM(:,2)= u_3211; 



% Noise

%noise.Enabler = 0;
noise.Enabler = 1;

noise.pos_stand_dev = noise.Enabler * 0.0011;                            	%[m]

noise.vel_stand_dev = noise.Enabler * 0.01;                               %[m/s]

noise.acc_stand_dev = noise.Enabler * 0.5;                               %[m/s^2]

noise.attitude_stand_dev = noise.Enabler * deg2rad(0.33);                 %[rad]
noise.ang_rate_stand_dev = noise.Enabler * deg2rad(1);                   %[rad/s]

% Delays

delay.position_filter = 1;
delay.attitude_filter = 1;
delay.mixer = 1;

% Load controller parameters

parameters_controller                    

% M injection example (sweeep: first column time vector, secondo column time history of pitching moment) 

SetPoint=[0,0];

% Values selected

t=ExcitationM(:,1);

simulation_time=t(end)-t(1);

theta_true = [0.1068 0.1192 -5.9755 -2.6478 -10.1647 450.71];
% Launch SIMULATOR TRUE VALUE
Xu = theta_true(1);
Xq = theta_true(2);
Mu = theta_true(3);
Mq = theta_true(4);
Xd = theta_true(5);
Md = theta_true(6);



A=[Xu, Xq, -9.81; Mu, Mq, 0; 0, 1, 0];

B=[Xd; Md; 0];

C=[1, 0, 0; 0, 1, 0; 0, 0, 1; Xu, Xq, 0]; 

D=[0; 0; 0; Xd];


out_true = sim("Simulator_Single_Axis"); 


% Launch SIMULATOR ESTIMATED VALUE

Xu = theta_est(1);
Xq = theta_est(2);
Mu = theta_est(3);
Mq = theta_est(4);
Xd = theta_est(5);
Md = theta_est(6);


clear A B C D
A=[Xu, Xq, -9.81; Mu, Mq, 0; 0, 1, 0];

B=[Xd; Md; 0];

C=[1, 0, 0; 0, 1, 0; 0, 0, 1; Xu, Xq, 0]; 

D=[0; 0; 0; Xd];


out_est = sim("Simulator_Single_Axis"); 

figure
plot(out_true.q.Time,out_true.q.Data)
hold on
plot(out_est.q.Time,out_est.q.Data)
legend('q true','q estimated')

figure
plot(out_true.ax.Time,out_true.ax.Data)
hold on
plot(out_est.ax.Time,out_est.ax.Data)
legend('ax true','ax estimated')