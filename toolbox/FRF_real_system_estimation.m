

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


%% STEP 6: VISUALIZATION (dB + deg)
fprintf('\nStep 6: Plotting results...\n');

% FIGURA 1: FRF Bode Plots
figure('Position', [100,100,1200,800], 'Name', 'FRF Welch - Bode');

subplot(2,2,1); 
semilogx(f, 20*log10(abs(H_q)), 'b-', 'LineWidth', 2); 
grid on; title('Pitch Rate |H_q| '); ylabel('[dB]'); 
ylim([-60 20]); xlim([f(1) f(end)]);

subplot(2,2,2); 
semilogx(f, 20*log10(abs(H_ax)), 'b-', 'LineWidth', 2);
grid on; title('Acceleration |H_{ax}| '); ylabel('[dB]'); 
ylim([-60 40]); xlim([f(1) f(end)]);

subplot(2,2,3); 
semilogx(f, angle(H_q)*180/pi, 'r-', 'LineWidth', 2);
grid on; xlabel('Frequency [Hz]'); ylabel('Phase [°]'); 
ylim([-180 180]); xlim([f(1) f(end)]);

subplot(2,2,4); 
semilogx(f, angle(H_ax)*180/pi, 'r-', 'LineWidth', 2);
grid on; xlabel('Frequency [Hz]'); ylabel('Phase [°]'); 
ylim([-180 180]); xlim([f(1) f(end)]);

sgtitle(sprintf('FRF Estimation - Welch Method )'), 'FontSize', 14, 'FontWeight', 'bold');

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




end

