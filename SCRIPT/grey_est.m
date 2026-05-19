


function [theta_est,std_theta,cov_theta] = grey_est (sigma_q,sigma_ax,f_min,f_max)

% Noise
%noise.Enabler = 0;
noise.Enabler = 1;
noise.pos_stand_dev = noise.Enabler * 0.0011;                            %[m]
noise.vel_stand_dev = noise.Enabler * 0.01;                               %[m/s]
noise.attitude_stand_dev = noise.Enabler * deg2rad(0.33);                 %[rad]
% noise su accellerometro
noise.acc_stand_dev = noise.Enabler * sigma_ax;                               %[m/s^2]
% noise su q
noise.ang_rate_stand_dev = noise.Enabler * deg2rad(sigma_q);                   %[rad/s]

% bisogna passare tutto noise non solo il campo modificato
assignin('base', 'noise',noise);
assignin('base', 'noise',noise);

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
    fprintf('(fs=%.1f Hz, %d samples per run)', fs, samples_per_run);

    % Fourier transform
    y = [q, ax];
    N = length(y);
    Y = fft(y);
    U = fft(u);

%% REAL SYSTEM FRF
% [H_q_real_sys,H_ax_real_sys,freq_vect] = FRF_real_system_estimation(u,q,ax,t,fs);
 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AUTOMATIC APPROACH %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% =============================================
%  STEP 0: Generiamo i dati
%  =============================================

% Procedure: Isolate dynamic response by removing the mean
    % Rationale: Spectral analysis focuses on fluctuations. Removing the DC
    % offset prevents a large, distorting spike at 0 Hz in the spectrum.
    u_no_mean  = u - mean(u);
    q_no_mean  = q - mean(q);
    ax_no_mean = ax - mean(ax);
%-->Step#3:Spectral Estimation (Welch's Method)
    % Apply the Welch method to compute smooth estimates of the auto- and
    % cross-power spectral densities (PSDs). This robustly handles the high
    % variance of a simple periodogram.
    fprintf('Step 3: Performing spectral estimation using Welch''s method... (σ_q=%.2f°/s, σ_ax=%.4f m/s²)\n', sigma_q, sigma_ax);
    % Define Welch method parameters
    M = 2048;               % Segment length (a power of 2 for FFT efficiency)
    noverlap = M / 2;       % 50% overlap between segments
    win = hanning(M);       % Hanning window to reduce spectral leakage
    % Compute the required spectra using MATLAB's built-in functions
    % 'pwelch' and 'cpsd' implement the entire 4-step Welch procedure:
    % (Segment -> Window -> Transform -> Average)
    % Input Autospectrum (Suu)
    [Suu, f] = pwelch(u_no_mean, win, noverlap, [], fs);
    % Cross-spectrum between input (u) and pitch rate (q)
    [Suq, ~] = cpsd(u_no_mean, q_no_mean, win, noverlap, [], fs);
    % Cross-spectrum between input (u) and longitudinal acceleration (ax)
    [Suax, ~] = cpsd(u_no_mean, ax_no_mean, win, noverlap, [], fs);
%  %-->Step#4:FRF Computation 
%  % Compute the non-parametric FRF estimate for each output channel using
%  % the fundamental spectral relationship: G(f) = S_uy(f) / S_uu(f).
% 
% 
% 
%-->Step#5:Quality Asesment   
    % Compute the coherence function for each input-output pair to
    % quantitatively assess the reliability of the FRF estimate.
    fprintf('Step 5: Assessing FRF quality with the coherence function... (σ_q=%.2f°/s, σ_ax=%.4f m/s²)\n', sigma_q, sigma_ax);
    % To compute coherence, we first need the output autospectra (Sqq, Saxax)
    [Sqq, ~] = pwelch(q_no_mean, win, noverlap, [], fs);
    [Saxax, ~] = pwelch(ax_no_mean, win, noverlap, [], fs);
    % Coherence formula: gamma^2 = |S_uy|^2 / (S_uu * S_yy)
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
    xline(f_min, 'r--', 'LineWidth', 2.5, 'Label', sprintf('f_{min}=%.2f Hz', f_min));
    xline(f_max, 'g--', 'LineWidth', 2.5, 'Label', sprintf('f_{max}=%.2f Hz', f_max));
    % Longitudinal Acceleration Channel
    subplot(2,2,2);
    semilogx(f, 20*log10(abs(Suax)));
    grid on; title(sprintf('CSPD ax (σ_q=%.2f°/s, σ_ax=%.4f m/s²)', sigma_q, sigma_ax));
    xline(f_min, 'r--', 'LineWidth', 2.5, 'Label', sprintf('f_{min}=%.2f Hz', f_min));
    xline(f_max, 'g--', 'LineWidth', 2.5, 'Label', sprintf('f_{max}=%.2f Hz', f_max));
    % Pitch Rate coherence
    subplot(2,2,3);
    semilogx(f, gamma_sq_q);
    hold on;
    plot([f(1), f(end)], [0.7, 0.7], 'r--', 'LineWidth', 1.5); % Threshold line
    grid on; title(sprintf('Coherence: Input (u) to Pitch Rate (q) (σ_q=%.2f°/s, σ_ax=%.4f m/s²)', sigma_q, sigma_ax));
    ylim([0 1.05]); ylabel('Coherence \gamma^2'); legend('Coherence', 'Trust Threshold');
    xline(f_min, 'r--', 'LineWidth', 2.5, 'Label', sprintf('f_{min}=%.2f Hz', f_min));
    xline(f_max, 'g--', 'LineWidth', 2.5, 'Label', sprintf('f_{max}=%.2f Hz', f_max));
    % Longitudinal Acceleration coherence
    subplot(2,2,4);
    semilogx(f, gamma_sq_ax);
    hold on;
    plot([f(1), f(end)], [0.7, 0.7], 'r--', 'LineWidth', 1.5); % Threshold line
    grid on; title(sprintf('Coherence: Input (u) to Longitudinal Accel. (ax) (σ_q=%.2f°/s, σ_ax=%.4f m/s²)', sigma_q, sigma_ax));
    ylim([0 1.05]); xlabel('Frequency [Hz]'); ylabel('Coherence \gamma^2');
    xline(f_min, 'r--', 'LineWidth', 2.5, 'Label', sprintf('f_{min}=%.2f Hz', f_min));
    xline(f_max, 'g--', 'LineWidth', 2.5, 'Label', sprintf('f_{max}=%.2f Hz', f_max));


%% =============================================
%  STEP 5: CREA MODELLO INIZIALE
%  =============================================
aux = {};
% Ordine dei parametri: [Xu; Xq; Mu; Mq; Xd; Md]
theta_0 = [-0.001; 0.001; -5.8; -4.7; -10; 120];

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
% Vettore di frequenze (in Hz, non rad/s per semplicità)
freq_Hz = (0:N-1)' / (N * Ts);

% =============================================
%  STEP 2: SELEZIONA INTERVALLO DI FREQUENZE
%  =============================================
% Definisci il range di interesse (in Hz)
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

data = iddata(Y_filtered, U_filtered,0.004, 'Frequency', W_filtered);
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
opt.SearchOptions.MaxIterations = 100;
% assegna un range in cui cercare i parametri
theta_true = [-0.1068 0.1192 -5.9755 -2.6478 -10.1647 450.71];
% for i = 1:length(theta_true)
%     if theta_true(i) < 0
%         init_sys.Structure.Parameters(1).Minimum(i) = theta_true(i)* 2;
%         init_sys.Structure.Parameters(1).Maximum(i) = theta_true(i) /2;
%     else
%         init_sys.Structure.Parameters(1).Minimum(i) = theta_true(i) /2;
%         init_sys.Structure.Parameters(1).Maximum(i) = theta_true(i) *2;
%     end
% end
%% =============================================
%  STEP 7: ESEGUI GREYEST CON DATI FREQUENZIALI
%  =============================================
fprintf('\n========== RUNNING GREYEST (FREQUENCY DOMAIN ESTIMATION) (σ_q=%.2f°/s, σ_ax=%.4f m/s²) ==========\n', sigma_q, sigma_ax);
[sys, x0] = greyest(data, init_sys, opt);
fprintf('\n========== ESTIMATION COMPLETED (σ_q=%.2f°/s, σ_ax=%.4f m/s²) ==========\n', sigma_q, sigma_ax);
%% =============================================
%  STEP 8: ANALIZZA RISULTATI
%  =============================================
fprintf('\n========== ESTIMATED PARAMETERS (σ_q=%.2f°/s, σ_ax=%.4f m/s²) ==========\n', sigma_q, sigma_ax);

% PARAMITER ESTIMATE
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
fprintf('\n========== FIT REPORT (σ_q=%.2f°/s, σ_ax=%.4f m/s²) ==========\n', sigma_q, sigma_ax);
disp(sys.Report.Fit);
fprintf('\n========== COVARIANCE & UNCERTAINTY (σ_q=%.2f°/s, σ_ax=%.4f m/s²) ==========\n', sigma_q, sigma_ax);
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
z_2sigma = 1.960;   % 95% CI 
z_3sigma = 2.576;   % 99% CI
fprintf('\n95%% Confidence Intervals (σ_q=%.2f°/s, σ_ax=%.4f m/s²):\n', sigma_q, sigma_ax);
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
% T1.2 -  Validation Against True Values
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
%   CONFRONTO TRANSFER FUNCTIONS: TRUE vs ESTIMATED
%  =============================================

% True model
theta_true_col = [-0.1068; 0.1192; -5.9755; -2.6478; -10.1647; 450.71];
true_sys = idgrey(@System_matrix, theta_true_col, 'c', {}, 0);

% Frequenze in rad/s
w = logspace(-2, 2, 500);

% Bode true model
[mag_true, phase_true] = bode(true_sys, w);

% Bode estimated model
[mag_est, phase_est] = bode(sys, w);

% Rimuove dimensioni singleton
mag_true   = squeeze(mag_true);
phase_true = squeeze(phase_true);

mag_est   = squeeze(mag_est);
phase_est = squeeze(phase_est);

% Per sistema 2 output - 1 input:
% mag(1,:) = output q
% mag(2,:) = output ax

figure('Position', [150, 150, 1000, 750], ...
       'Name', 'Transfer Function Comparison: Real vs Estimated');

sgtitle('Transfer Function Comparison: Real vs Estimated', ...
        'FontSize', 14, 'FontWeight', 'bold');

% =========================
% Pitch rate - Magnitude
% =========================
subplot(2,2,1)
semilogx(w, 20*log10(abs(mag_true(1,:))), 'b-', 'LineWidth', 1.5);
hold on
semilogx(w, 20*log10(abs(mag_est(1,:))), 'r--', 'LineWidth', 1.5);
grid on
ylabel('Magnitude [dB]')
title('Pitch rate q')
legend('True Model', 'Identified Model', 'Location', 'best')

% =========================
% Longitudinal acceleration - Magnitude
% =========================
subplot(2,2,2)
semilogx(w, 20*log10(abs(mag_true(2,:))), 'b-', 'LineWidth', 1.5);
hold on
semilogx(w, 20*log10(abs(mag_est(2,:))), 'r--', 'LineWidth', 1.5);
grid on
ylabel('Magnitude [dB]')
title('Longitudinal acceleration a_x')
legend('True Model', 'Identified Model', 'Location', 'best')

% =========================
% Pitch rate - Phase
% =========================
subplot(2,2,3)
semilogx(w, phase_true(1,:), 'b-', 'LineWidth', 1.5);
hold on
semilogx(w, phase_est(1,:), 'r--', 'LineWidth', 1.5);
grid on
xlabel('Frequency [rad/s]')
ylabel('Phase [deg]')
title('Pitch rate q')

% =========================
% Longitudinal acceleration - Phase
% =========================
subplot(2,2,4)
semilogx(w, phase_true(2,:), 'b-', 'LineWidth', 1.5);
hold on
semilogx(w, phase_est(2,:), 'r--', 'LineWidth', 1.5);
grid on
xlabel('Frequency [rad/s]')
ylabel('Phase [deg]')
title('Longitudinal acceleration a_x')
% %%
% % =============================================
% %   CONFRONTO MODELLO vs DATI FREQUENZIALI
% %  =============================================
% figure('Position', [100, 100, 1400, 900], 'Name', sprintf('Compare: Data vs Model (σ_q=%.2f°/s, σ_ax=%.4f m/s²)', sigma_q, sigma_ax));
% 
% [mag, phase, wout] = bode(sys);
% grid on;
% title(sprintf('Bode Plot of Identified Grey-Box Model (σ_q=%.2f°/s, σ_ax=%.4f m/s²)', sigma_q, sigma_ax), 'FontSize', 14);
% 
% opt_comp = compareOptions('InitialCondition', 'zero');
% compare(data, sys, opt_comp);
% % Post-processing per rendere leggibile
% h = gcf;
% for i = 1:length(h.Children)
%     ax = h.Children(i);
%     if isa(ax, 'matlab.graphics.axis.Axes')
%         ax.FontSize = 10;
%         ax.YLabel.FontSize = 10;
%         ax.XLabel.FontSize = 10;
%         ax.Title.FontSize = 11;
%     end
% end
% 
% % Estrazione dati dal bode (reshape da 3D a 2D)
% mag_q  = squeeze(mag(1, 1, :));    % Pitch rate magnitude
% mag_ax = squeeze(mag(2, 1, :));    % Acceleration magnitude
% phase_q  = squeeze(phase(1, 1, :)); % Pitch rate phase
% phase_ax = squeeze(phase(2, 1, :)); % Acceleration phase
% freq_Hz = wout / (2*pi);            % Converti da rad/s a Hz
% 
% figure('Position', [100, 100, 1200, 800], 'Name', sprintf('FRF Welch - Bode (σ_q=%.2f°/s, σ_ax=%.4f m/s²)', sigma_q, sigma_ax));
% 
% % Magnitude plots
% subplot(2, 2, 1); 
% semilogx(freq_vect, 20*log10(abs(H_q_real_sys)), 'b-', 'LineWidth', 2); 
% hold on
% semilogx(freq_Hz, 20*log10(mag_q), 'r-', 'LineWidth', 2); 
% grid on; 
% title(sprintf('Pitch Rate |H_q| (σ_q=%.2f°/s, σ_ax=%.4f m/s²)', sigma_q, sigma_ax)); 
% ylabel('[dB]'); 
% xlim([10^-1 10^2]);
% legend('real system', 'estimated')
% 
% subplot(2, 2, 2); 
% semilogx(freq_vect, 20*log10(abs(H_ax_real_sys)), 'b-', 'LineWidth', 2);
% hold on
% semilogx(freq_Hz, 20*log10(mag_ax), 'r-', 'LineWidth', 2); 
% grid on; 
% title(sprintf('Acceleration |H_{ax}| (σ_q=%.2f°/s, σ_ax=%.4f m/s²)', sigma_q, sigma_ax)); 
% ylabel('[dB]'); 
% xlim([10^-1 10^2]);
% legend('real system', 'estimated')
% 
% % Phase plots
% subplot(2, 2, 3); 
% semilogx(freq_vect, angle(H_q_real_sys)*180/pi, 'b-', 'LineWidth', 2);
% hold on
% semilogx(freq_Hz, phase_q, 'r-', 'LineWidth', 2);
% grid on; 
% xlabel('Frequency [Hz]'); 
% ylabel('Phase [°]'); 
% title(sprintf('Pitch Rate Phase q (σ_q=%.2f°/s, σ_ax=%.4f m/s²)', sigma_q, sigma_ax))
% xlim([10^-1 10^2]);
% legend('real system', 'estimated')
% 
% subplot(2, 2, 4); 
% semilogx(freq_vect, angle(H_ax_real_sys)*180/pi, 'b-', 'LineWidth', 2);
% hold on
% semilogx(freq_Hz, phase_ax, 'r-', 'LineWidth', 2);
% grid on; 
% xlabel('Frequency [Hz]'); 
% ylabel('Phase [°]');
% title(sprintf('Acceleration Phase ax (σ_q=%.2f°/s, σ_ax=%.4f m/s²)', sigma_q, sigma_ax))
% xlim([10^-1 10^2]);
% legend('real system', 'estimated')
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


% Plot input sequence
figure('Name', sprintf('3211 Sequence Input (sigma_q=%.2f deg/s, sigma_ax=%.4f m/s^2)', ...
       sigma_q, sigma_ax), ...
       'Position', [200, 200, 900, 400]);

stairs(t_3211, u_3211, 'b-', 'LineWidth', 2.0);
grid on;

xlabel('Time [s]', 'FontSize', 12);
ylabel('Input amplitude', 'FontSize', 12);
title('3211 Sequence (\tau = 2 s, A = 0.1)', ...
      'FontSize', 14, 'FontWeight', 'bold');

xlim([0 Tend]);
ylim([-0.12 0.12]);

set(gca, 'FontSize', 11);


% % (facoltativo) visualizza
% figure('Name', sprintf('3211 Sequence Input (σ_q=%.2f°/s, σ_ax=%.4f m/s²)', sigma_q, sigma_ax)); 
% stairs(t_3211, u_3211, 'LineWidth', 2); 
% grid on;
% xlabel('Time [s]'); 
% ylabel('Input amplitude'); 
% title(sprintf('3211 Sequence (σ_q=%.2f°/s, σ_ax=%.4f m/s²)', sigma_q, sigma_ax));

ExcitationM(:,1)= t_3211; 
ExcitationM(:,2)= u_3211; 
t=ExcitationM(:,1);
simulation_time=t(end)-t(1);

% il modello simulink prende i dati dal base workspace
evalin('base', 'clear ExcitationM');
assignin('base', 'ExcitationM', ExcitationM); % sovrascrive la variabile nel base workspace con quella della funzione
assignin('base', 'simulation_time',simulation_time);

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

assignin('base', 'A',A);
assignin('base', 'B',B);
assignin('base', 'C',C);
assignin('base', 'D',D);
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

assignin('base', 'A',A);
assignin('base', 'B',B);
assignin('base', 'C',C);
assignin('base', 'D',D);

out_est = sim("Simulator_Single_Axis"); 

% PLOT COME NELL'IMMAGINE - Layout 2x2
figure('Name', sprintf('Figure 1: 3211 sequence validation (σ_q=%.2f°/s, σ_ax=%.4f m/s²)', sigma_q, sigma_ax));
% Subplot 1: q response
subplot(2,2,1)
plot(out_true.q.Time, out_true.q.Data, 'b-', 'LineWidth', 1.5)
hold on
plot(out_est.q.Time, out_est.q.Data, 'r-', 'LineWidth', 1.5)
grid on
xlabel('time')
ylabel('q')
title(sprintf('q response to 3211 sequence (σ_q=%.2f°/s, σ_ax=%.4f m/s²)', sigma_q, sigma_ax))
legend('Real', 'Estimated', 'Location', 'best')

% Subplot 2: ax response
subplot(2,2,2)
plot(out_true.ax.Time, out_true.ax.Data, 'b-', 'LineWidth', 1.5)
hold on
plot(out_est.ax.Time, out_est.ax.Data, 'r-', 'LineWidth', 1.5)
grid on
xlabel('time')
ylabel('ax')
title(sprintf('ax response to 3211 sequence (σ_q=%.2f°/s, σ_ax=%.4f m/s²)', sigma_q, sigma_ax))
legend('Real', 'Estimated', 'Location', 'best')

% Subplot 3: q error
subplot(2,2,3)
q_error = out_true.q.Data - out_est.q.Data;
plot(out_true.q.Time, q_error, 'b-', 'LineWidth', 1.5)
grid on
xlabel('time')
ylabel('q error')
title(sprintf('q error (σ_q=%.2f°/s, σ_ax=%.4f m/s²)', sigma_q, sigma_ax))

% Subplot 4: ax error
subplot(2,2,4)
ax_error = out_true.ax.Data - out_est.ax.Data;
plot(out_true.ax.Time, ax_error, 'b-', 'LineWidth', 1.5)
grid on
xlabel('time')
ylabel('ax error')
title(sprintf('ax error (σ_q=%.2f°/s, σ_ax=%.4f m/s²)', sigma_q, sigma_ax))

%% CALCOLO DEL FITTING NEL TEMPO 

% Recupero le varianze vere
var_q = noise.ang_rate_stand_dev^2; 
var_ax = noise.acc_stand_dev^2;

% Errore Pitch Rate (q)
e_q = out_true.q.Data - out_est.q.Data;
% Calcolo J 
J_q = sum(e_q.^2) / (2 * var_q); 
% Calcolo Fit Percentuale
fit_q_percent = 100 * (1 - norm(e_q) / norm(out_true.q.Data - mean(out_true.q.Data)));

% 2. Errore Accelerazione (ax)
e_ax = out_true.ax.Data - out_est.ax.Data;
% Calcolo J 
J_ax = sum(e_ax.^2) / (2 * var_ax);
% Calcolo Fit Percentuale
fit_ax_percent = 100 * (1 - norm(e_ax) / norm(out_true.ax.Data - mean(out_true.ax.Data)));

% 3. Risultati a schermo
fprintf('\n========== VALIDAZIONE NEL TEMPO (Sequenza 3211) ==========\n');
fprintf('Pitch Rate (q) : J (OE) = %10.2f  |  Fit = %.2f%%\n', J_q, fit_q_percent);
fprintf('Accel. (ax)    : J (OE) = %10.2f  |  Fit = %.2f%%\n', J_ax, fit_ax_percent);
fprintf('===========================================================\n');

end
