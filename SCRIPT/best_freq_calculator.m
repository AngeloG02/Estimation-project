function [f_min,f_max] = best_freq_calculator (sigma_q,sigma_ax)

% Noise
%noise.Enabler = 0;
noise.Enabler = 1;
noise.pos_stand_dev = noise.Enabler * 0.0011;                            %[m]
noise.vel_stand_dev = noise.Enabler * 0.01;                              %[m/s]
noise.attitude_stand_dev = noise.Enabler * deg2rad(0.33);                %[rad]
% noise su accelerometro
noise.acc_stand_dev = noise.Enabler * sigma_ax;                          %[m/s^2]
% noise su q
noise.ang_rate_stand_dev = noise.Enabler * deg2rad(sigma_q);             %[rad/s]

% bisogna passare tutto noise non solo il campo modificato
assignin('base', 'noise', noise);

%% STEP 1: DATA ACQUISITION
fprintf('\nStep 1: Acquiring data from Simulink runs... (σ_q=%.2f°/s, σ_ax=%.4f m/s²)\n', sigma_q, sigma_ax);

% Run simulation
out = sim("Simulator_Single_Axis");

% Extract signals from this Simulink run
t = out.q.Time;
u = out.Mtot.Data;      % Input: Momento totale
q = out.q.Data;         % Output 1: Pitch rate
ax = out.ax.Data;       % Output 2: Accel. longitudinale

Ts = t(2) - t(1);
fs = 1/Ts;
samples_per_run = length(t);
fprintf('(fs=%.1f Hz, %d samples per run)\n', fs, samples_per_run);

% Fourier transform - SOLO FREQUENZE POSITIVE
y = [q, ax];

N = length(u);
Y = fft(y);
U = fft(u);

Npos = floor(N/2) + 1;
Y = Y(1:Npos, :);
U = U(1:Npos);

freq_Hz = (0:Npos-1)' * fs / N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MANUAL APPROACH %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MODEL SYSTEM FRF
% G = @(theta,freq_vect)Model_system_FRF(theta,freq_vect);

%% REAL SYSTEM FRF
[H_q_real_sys,H_ax_real_sys,freq_vect] = FRF_real_system_estimation(u,q,ax,t,fs);
% G_m = [H_q_real_sys,H_ax_real_sys];

%% OTTIMIZZATION
% theta_0 = [-0.5; -0.1; 0.5; -20; -1; 1];
% OPTIMIZATION_WITH_FMINCON
% [theta_hat_oe, J_hat_oe] = outputError_FD(G_m, freq_vect, G, theta_0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AUTOMATIC APPROACH %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% =============================================
%  STEP 0: Generiamo i dati
% =============================================

% Procedure: Isolate dynamic response by removing the mean
u_no_mean  = u - mean(u);
q_no_mean  = q - mean(q);
ax_no_mean = ax - mean(ax);

%-->Step#3: Spectral Estimation (Welch's Method)
fprintf('Step 3: Performing spectral estimation using Welch''s method... (σ_q=%.2f°/s, σ_ax=%.4f m/s²)\n', sigma_q, sigma_ax);

M = 2048;               % Segment length
noverlap = M / 2;       % 50% overlap
win = hanning(M);       % Hanning window

% Input Autospectrum (Suu)
[Suu, f] = pwelch(u_no_mean, win, noverlap, [], fs);

% Cross-spectrum between input (u) and pitch rate (q)
[Suq, ~] = cpsd(u_no_mean, q_no_mean, win, noverlap, [], fs);

% Cross-spectrum between input (u) and longitudinal acceleration (ax)
[Suax, ~] = cpsd(u_no_mean, ax_no_mean, win, noverlap, [], fs);

%-->Step#5: Quality Assessment
fprintf('Step 5: Assessing FRF quality with the coherence function... (σ_q=%.2f°/s, σ_ax=%.4f m/s²)\n', sigma_q, sigma_ax);

[Sqq, ~] = pwelch(q_no_mean, win, noverlap, [], fs);
[Saxax, ~] = pwelch(ax_no_mean, win, noverlap, [], fs);

gamma_sq_q = (abs(Suq).^2) ./ (Suu .* Sqq);
gamma_sq_ax = (abs(Suax).^2) ./ (Suu .* Saxax);

%-->Step#6:Vizualization   
    % Plot the results for analysis.
    fprintf('Step 6: Generating plots... (σ_q=%.2f°/s, σ_ax=%.4f m/s²)\n', sigma_q, sigma_ax);
    % Plot 1: Estimated FRFs (Bode Plot)
    figure('Name', sprintf('Spectral Analysis (σ_q=%.2f°/s, σ_ax=%.4f m/s²)', sigma_q, sigma_ax))
    % Pitch Rate Channel
    subplot(2,2,1);
    semilogx(f, 20*log10(abs(Suq)));
    grid on; title(sprintf('CSPD q (σ_q=%.2f°/s, σ_ax=%.4f m/s²)', sigma_q, sigma_ax)); ylabel('Magnitude [dB]');
    % Longitudinal Acceleration Channel
    subplot(2,2,2);
    semilogx(f, 20*log10(abs(Suax)));
    grid on; title(sprintf('CSPD ax (σ_q=%.2f°/s, σ_ax=%.4f m/s²)', sigma_q, sigma_ax));
    % Pitch Rate coherence
    subplot(2,2,3);
    semilogx(f, gamma_sq_q);
    hold on;
    plot([f(1), f(end)], [0.7, 0.7], 'r--', 'LineWidth', 1.5); % Threshold line
    grid on; title(sprintf('Coherence: Input (u) to Pitch Rate (q) (σ_q=%.2f°/s, σ_ax=%.4f m/s²)', sigma_q, sigma_ax));
    ylim([0 1.05]); ylabel('Coherence \gamma^2'); legend('Coherence', 'Trust Threshold');
    % Longitudinal Acceleration coherence
    subplot(2,2,4);
    semilogx(f, gamma_sq_ax);
    hold on;
    plot([f(1), f(end)], [0.7, 0.7], 'r--', 'LineWidth', 1.5); % Threshold line
    grid on; title(sprintf('Coherence: Input (u) to Longitudinal Accel. (ax) (σ_q=%.2f°/s, σ_ax=%.4f m/s²)', sigma_q, sigma_ax));
    ylim([0 1.05]); xlabel('Frequency [Hz]'); ylabel('Coherence \gamma^2');



%% =============================================
%  STEP 5: CREA MODELLO INIZIALE
% =============================================
aux = {};
theta_0 = [0.1; 0.1; -0.1; -0.1; -1;100 ];
T = 0;

init_sys = idgrey(@System_matrix, theta_0, 'c', aux, T);

% Controlla stabilità
eig_init = eig(init_sys.A);
fprintf('\n========== INITIAL MODEL CHECK (σ_q=%.2f°/s, σ_ax=%.4f m/s²) ==========\n', sigma_q, sigma_ax);
fprintf('Eigenvalues of A:\n');
disp(eig_init);

if all(real(eig_init) < 0)
    fprintf('✓ Model is STABLE\n');
else
    fprintf('✗ WARNING: Model is UNSTABLE!\n');
    fprintf('  Adjust theta_0 values.\n');
end

%% =============================================
%  STEP 6: CONFIGURA OPZIONI
% =============================================
opt = greyestOptions;
opt.InitialState = 'backcast';
opt.Focus = 'simulation';
opt.Display = 'on';
opt.EnforceStability = false;
opt.SearchMethod = 'auto';
opt.SearchOptions.MaxIterations = 100;

% assegna un range in cui cercare i parametri
theta_true = [-0.1068 0.1192 -5.9755 -2.6478 -10.1647 450.71];
% for i = 1:length(theta_true)
%     if theta_true(i) < 0
%         init_sys.Structure.Parameters(1).Minimum(i) = theta_true(i) * 2;
%         init_sys.Structure.Parameters(1).Maximum(i) = theta_true(i) / 2;
%     else
%         init_sys.Structure.Parameters(1).Minimum(i) = theta_true(i) / 2;
%         init_sys.Structure.Parameters(1).Maximum(i) = theta_true(i) * 2;
%     end
% end

%% =============================================
%  STEP 6.5: TROVA BANDA OTTIMA AUTOMATICAMENTE
% =============================================
result_band = optimize_fband(u, q, ax, Ts, gamma_sq_q, gamma_sq_ax, f, init_sys, opt);

f_min = result_band.f_min;
f_max = result_band.f_max;

fprintf('\nOptimal frequency range found: %.2f - %.2f Hz\n', f_min, f_max);
fprintf('Best internal fit_q  = %.2f %%\n', result_band.fit_q);
fprintf('Best internal fit_ax = %.2f %%\n', result_band.fit_ax);
fprintf('Best internal fit_sum = %.2f %%\n', result_band.fit_sum);



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

end