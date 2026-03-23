

function [H_q,H_ax,f] = FRF_real_system_estimation(u_combined,q_combined,ax_combined,t_combined,fs)
set(0,'DefaultFigureWindowStyle','docked')


u_raw = u_combined;
q_raw = q_combined;
ax_raw = ax_combined;
N_total = length(u_raw);
t = t_combined;

%% STEP 2: PRE-PROCESSING (Zero-mean)
fprintf('\nStep 2: Zero-mean correction...\n');
u = u_raw - mean(u_raw);
q = q_raw - mean(q_raw);
ax = ax_raw - mean(ax_raw);
fprintf('  ✓ Signals zero-meaned (DC removed)\n');

%% STEP 3: SPECTRAL ESTIMATION (Welch's Method)
fprintf('\nStep 3: Welch spectral estimation...\n');

% Parametri Welch ottimizzati
M = 2048;           % Lunghezza segmento (power of 2)
noverlap = M/2;     % 50% overlap
win = hanning(M);  % BARTLETT WINDOW (good sidelobe suppression)

% Number of segments available
K = floor((N_total - M) / (M - noverlap)) + 1;

%% STEP 4: FRF COMPUTATION (tfestimate - H1 Estimator)
fprintf('\nStep 4: Computing FRF with tfestimate...\n');

% ---> AGGIUNGIAMO NFFT PER AVERE PIÙ PUNTI <---
% Usa potenze di 2. 16384 o 32768 ti daranno una griglia fittissima.
NFFT = 16384; 

% H1 Estimator: assumes noise uncorrelated with input (most common)
% Inseriamo NFFT al posto di "[]"
[H_q, f] = tfestimate(u, q, win, noverlap, NFFT, fs, 'Estimator', 'H1');
[H_ax, ~] = tfestimate(u, ax, win, noverlap, NFFT, fs, 'Estimator', 'H1');
fprintf('\nStep 4: Computing FRF with tfestimate...\n');


% % H1 Estimator: assumes noise uncorrelated with input (most common)
% [H_q, f] = tfestimate(u, q, win, noverlap, [], fs, 'Estimator', 'H1');
% [H_ax, ~] = tfestimate(u, ax, win, noverlap, [], fs, 'Estimator', 'H1');
% fprintf('\nStep 4: Computing FRF with tfestimate...\n');
% 

% %% STEP 5: QUALITY ASSESSMENT (mscohere)
% fprintf('\nStep 5: Coherence analysis with mscohere...\n');
% 
% [coh2_q, ~] = mscohere(u, q, win, noverlap, [], fs);
% [coh2_ax, ~] = mscohere(u, ax, win, noverlap, [], fs);
% 
% coh_q_mean = mean(coh2_q);
% coh_ax_mean = mean(coh2_ax);

%% da scommentare

% %% STEP 6: VISUALIZATION (dB + deg)
% fprintf('\nStep 6: Plotting results...\n');
% 
% % FIGURA 1: FRF Bode Plots
% figure('Position', [100,100,1200,800], 'Name', 'FRF Welch - Bode');
% 
% subplot(2,2,1); 
% semilogx(f, 20*log10(abs(H_q)), 'b-', 'LineWidth', 2); 
% grid on; title('Pitch Rate |H_q| '); ylabel('[dB]'); 
% ylim([-60 20]); xlim([f(1) f(end)]);
% 
% subplot(2,2,2); 
% semilogx(f, 20*log10(abs(H_ax)), 'b-', 'LineWidth', 2);
% grid on; title('Acceleration |H_{ax}| '); ylabel('[dB]'); 
% ylim([-60 40]); xlim([f(1) f(end)]);
% 
% subplot(2,2,3); 
% semilogx(f, angle(H_q)*180/pi, 'r-', 'LineWidth', 2);
% grid on; xlabel('Frequency [Hz]'); ylabel('Phase [°]'); 
% ylim([-180 180]); xlim([f(1) f(end)]);
% 
% subplot(2,2,4); 
% semilogx(f, angle(H_ax)*180/pi, 'r-', 'LineWidth', 2);
% grid on; xlabel('Frequency [Hz]'); ylabel('Phase [°]'); 
% ylim([-180 180]); xlim([f(1) f(end)]);
% 
% sgtitle(sprintf('FRF Estimation - Welch Method )'), 'FontSize', 14, 'FontWeight', 'bold');
% 
% % FIGURA 2: Coherence
% figure('Position', [200,200,1000,500], 'Name', 'Coherence Functions');
% 
% subplot(1,2,1); 
% semilogx(f, coh2_q, 'g-', 'LineWidth', 2);
% hold on; 
% yline(0.6, 'r--', 'LineWidth', 2, 'Label', 'γ²=0.6 (Trust Threshold)');
% grid on; title('Pitch Rate Coherence'); ylabel('γ²'); ylim([0 1.1]);
% xlim([f(1) f(end)]);
% legend('γ²', 'Trust Threshold', 'Mean', 'Location', 'best');
% 
% subplot(1,2,2); 
% semilogx(f, coh2_ax, 'g-', 'LineWidth', 2);
% hold on; 
% yline(0.6, 'r--', 'LineWidth', 2, 'Label', 'γ²=0.6 (Trust Threshold)');
% grid on; title('Acceleration Coherence'); xlabel('Frequency [Hz]'); 
% ylabel('γ²'); ylim([0 1.1]); xlim([f(1) f(end)]);
% legend('γ²', 'Trust Threshold', 'Mean', 'Location', 'best');


%% STEP 7: CONFIDENCE-BASED MASKING


% fprintf('\nStep 7: Quality filtering...\n');
% 
% trust_threshold = 0.5;  % Standard engineering practice
% 
% % Definisci la maschera logica combinata
% index = (coh2_q > trust_threshold) & (coh2_ax > trust_threshold);
% 
% % Applica la maschera ai vettori
% H_q_trusted  = H_q(index);
% H_ax_trusted = H_ax(index);
% freq_vect    = f(index);


end

% %% ========================================================================
% % T1.1 FRF ESTIMATION (WELCH METHOD) - LONGITUDINAL DYNAMICS
% % Versione COMATTA con pwelch/cpsd + Magnitude dB + Phase deg
% % ========================================================================
% clear all; close all; clc;
% set(0,'DefaultFigureWindowStyle','docked')
% 
% fprintf('\n%s\n', repmat('=',60,1));
% fprintf('FRF ESTIMATION (WELCH) - pwelch/cpsd VERSION\n');
% fprintf('%s\n\n', repmat('=',60,1));
% 
% %% STEP 1: DATA ACQUISITION
% fprintf('Step 1: Acquiring data...\n');
% main_quad_ANTX;
% t = out.q.Time;
% Ts = t(2) - t(1);
% fs = 1/Ts;
% 
% % Extract signals
% u_raw = out.Mtot.Data;      % Input: Momento totale
% q_raw = out.q.Data;         % Output 1: Pitch rate
% ax_raw = out.ax.Data;       % Output 2: Accel. longitudinale
% 
% N = length(t);
% fprintf('  Data: %d samples (%.2f s) @ %.1f Hz\n', N, t(end), fs);
% 
% %% STEP 2: PRE-PROCESSING (Zero-mean)
% fprintf('Step 2: Zero-mean correction...\n');
% u = u_raw - mean(u_raw);
% q = q_raw - mean(q_raw);
% ax = ax_raw - mean(ax_raw);
% 
% %% STEP 3: SPECTRAL ESTIMATION (Welch's Method)
% fprintf('Step 3: Welch spectral estimation...\n');
% % Parametri Welch ottimizzati
% M = 2048;           % Lunghezza segmento (power of 2)
% noverlap = M/2;     % 50% overlap
% % win = hanning(M);   % Hanning window
% win = bartlett(M);  %  BARTLETT WINDOW
% 
% % INPUT Autospectrum
% [Suu, f] = pwelch(u, win, noverlap, [], fs);
% 
% % CROSS-SPECTRA
% [Suq, ~] = cpsd(u, q, win, noverlap, [], fs);    % Input → Pitch rate
% [Suax, ~] = cpsd(u, ax, win, noverlap, [], fs);  % Input → Acceleration
% 
% % OUTPUT Autospectra (per coerenza)
% [Sqq, ~] = pwelch(q, win, noverlap, [], fs);
% [Saxax, ~] = pwelch(ax, win, noverlap, [], fs);
% 
% fprintf('  Spectra computed @ %d freq points\n', length(f));
% 
% %% STEP 4: FRF COMPUTATION
% fprintf('Step 4: Computing FRF...\n');
% H_q = Suq ./ Suu;       % FRF Pitch Rate
% H_ax = Suax ./ Suu;     % FRF Acceleration
% 
% %% STEP 5: QUALITY ASSESSMENT (Coherence)
% fprintf('Step 5: Coherence analysis...\n');
% coh2_q = abs(Suq).^2 ./ (Suu .* Sqq);       % γ² per pitch rate
% coh2_ax = abs(Suax).^2 ./ (Suu .* Saxax);   % γ² per acceleration
% 
% fprintf('  Coherence ranges:\n');
% fprintf('    Pitch:   %.3f ≤ γ² ≤ %.3f\n', min(coh2_q), max(coh2_q));
% fprintf('    Accel:   %.3f ≤ γ² ≤ %.3f\n', min(coh2_ax), max(coh2_ax));
% 
% %% STEP 6: VISUALIZATION (dB + deg)
% fprintf('Step 6: Plotting results...\n');
% 
% % FIGURA 1: FRF Bode Plots
% figure('Position', [100,100,1200,800], 'Name', 'FRF Welch - Bode');
% subplot(2,2,1); semilogx(f, 20*log10(abs(H_q)), 'b-', 'LineWidth', 2); 
% grid on; title('Pitch Rate |H_q|'); ylabel('[dB]'); ylim([-60 20]);
% 
% subplot(2,2,2); semilogx(f, 20*log10(abs(H_ax)), 'b-', 'LineWidth', 2);
% grid on; title('Acceleration |H_{ax}|'); ylabel('[dB]'); ylim([-60 40]);
% 
% subplot(2,2,3); semilogx(f, angle(H_q)*180/pi, 'r-', 'LineWidth', 2);
% grid on; xlabel('Frequency [Hz]'); ylabel('Phase [°]'); ylim([-180 180]);
% 
% subplot(2,2,4); semilogx(f, angle(H_ax)*180/pi, 'r-', 'LineWidth', 2);
% grid on; xlabel('Frequency [Hz]'); ylabel('Phase [°]'); ylim([-180 180]);
% 
% sgtitle('T1.1 FRF Estimation - Welch Method (pwelch/cpsd)', 'FontSize', 14, 'FontWeight', 'bold');
% 
% % FIGURA 2: Coherence
% figure('Position', [200,200,1000,500], 'Name', 'Coherence Functions');
% subplot(1,2,1); semilogx(f, coh2_q, 'g-', 'LineWidth', 2);
% hold on; yline(0.6, 'r--', 'LineWidth', 2, 'Label', 'γ²=0.6');
% grid on; title('Pitch Rate Coherence'); ylabel('γ²'); ylim([0 1.1]);
% legend('γ²', 'Trust Threshold');
% 
% subplot(1,2,2); semilogx(f, coh2_ax, 'g-', 'LineWidth', 2);
% hold on; yline(0.6, 'r--', 'LineWidth', 2, 'Label', 'γ²=0.6');
% grid on; title('Acceleration Coherence'); xlabel('Frequency [Hz]'); ylim([0 1.1]);
% 
% fprintf('\n%s\n', repmat('=',60,1));
% fprintf('✅ ANALYSIS COMPLETE - Check the figures!\n');
% fprintf('%s\n', repmat('=',60,1));


