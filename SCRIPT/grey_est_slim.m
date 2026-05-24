% function [theta_est,std_theta] = grey_est_slim (sigma_q,sigma_ax)

% clc; clear; close all

main_quad_ANTX % inizializza i parametri che servono al modello

%% ========================================================================
%  INITIALIZATION
% =========================================================================

f_min = f_min_dx;
f_max = f_max_dx;

% CAMBIARE IN BASE ALLA CONFIGURAZIONE MIGLIORE
sigma_q  = sigma_q_dx;     % [deg/s]
sigma_ax = sigma_ax_dx;    % [m/s^2]

% Noise
noise.Enabler = 1;
noise.pos_stand_dev      = noise.Enabler * 0.0011;
noise.vel_stand_dev      = noise.Enabler * 0.01;
noise.attitude_stand_dev = noise.Enabler * deg2rad(0.33);
noise.acc_stand_dev      = noise.Enabler * sigma_ax;
noise.ang_rate_stand_dev = noise.Enabler * deg2rad(sigma_q);

assignin('base', 'noise', noise);

T = 0;
aux = {};

% Parametri iniziali
theta_0 = [-0.001; 0.001; -5.8; -4.7; -10; 120];

init_sys = idgrey(@System_matrix, theta_0, 'c', aux, T);

% Opzioni greyest
opt = greyestOptions;
opt.InitialState = 'backcast';
opt.Focus = 'simulation';
opt.Display = 'off';
opt.EnforceStability = false;
opt.SearchMethod = 'auto';
opt.SearchOptions.MaxIterations = 100;

% True parameters: [Xu; Xq; Mu; Mq; Xd; Md]
theta_true = [-0.1068; 0.1192; -5.9755; -2.6478; -10.1647; 450.71];

% True model
true_sys = idgrey(@System_matrix, theta_true, 'c', aux, T);
true_ss  = ss(true_sys);

Gq_true  = true_ss(1,1);   % M -> q
Gax_true = true_ss(2,1);   % M -> ax

poles_true = pole(true_ss);

zeros_q_true  = get_siso_zeros(Gq_true);
zeros_ax_true = get_siso_zeros(Gax_true);

%% ========================================================================
%  STEP 1: DATA ACQUISITION - MONTE CARLO
% =========================================================================

N_realizations = 300;

theta_est = zeros(6, N_realizations);
std_theta = zeros(6, N_realizations);

theta_mc = zeros(6, N_realizations);

% Poli e zeri separati per i due canali SISO
poles_q_mc  = cell(1, N_realizations);
poles_ax_mc = cell(1, N_realizations);

zeros_q_mc  = cell(1, N_realizations);
zeros_ax_mc = cell(1, N_realizations);

% Modelli salvati per FRF Monte Carlo
sys_mc = cell(1, N_realizations);

%% CREA ARRAY DI SimulationInput

simIn(1:N_realizations) = Simulink.SimulationInput('Simulator_Single_Axis');

for i = 1:N_realizations

    simIn(i) = simIn(i).setVariable('Seed_pos',   i);
    simIn(i) = simIn(i).setVariable('Seed_vel',   i + N_realizations);
    simIn(i) = simIn(i).setVariable('Seed_theta', i + 2*N_realizations);
    simIn(i) = simIn(i).setVariable('Seed_q',     i + 3*N_realizations);
    simIn(i) = simIn(i).setVariable('Seed_ax',    i + 4*N_realizations);

end

%% ESEGUI SIMULAZIONI IN PARALLELO

simOut = parsim(simIn, ...
    'ShowProgress', 'on', ...
    'TransferBaseWorkspaceVariables', 'on');

%% PROCESSA I RISULTATI

parfor realization = 1:N_realizations

    fprintf('  Realization %2d/%d: Processing... ', realization, N_realizations);

    out = simOut(realization);

    t  = out.q.Time;
    u  = out.Mtot.Data;
    q  = out.q.Data;
    ax = out.ax.Data;

    % Rimozione media
    u  = u  - mean(u);
    q  = q  - mean(q);
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

    % Dati in frequency domain
    data = iddata(Y_filtered, U_filtered, Ts, 'Frequency', W_filtered);

    data.InputName = 'Pitch Moment';
    data.InputUnit = '-';

    data.OutputName = {'Pitch rate', 'Longitudinal acceleration'};
    data.OutputUnit = {'rad/s', 'm/s^2'};

    data.TimeUnit = 's';

    % Grey-box estimation
    [sys, ~] = greyest(data, init_sys, opt);

    % Parametri stimati
    theta_est(:, realization) = getpvec(sys);
    theta_mc(:, realization)  = getpvec(sys);

    % Covarianza parametri
    cov_theta = getcov(sys);
    std_theta(:, realization) = sqrt(diag(cov_theta));

    % Conversione in state-space standard
    sys_ss = ss(sys);

    % Transfer functions SISO separate
    G_q  = sys_ss(1,1);   % M -> q
    G_ax = sys_ss(2,1);   % M -> ax

    % Poli e zeri del canale M -> q
    poles_q_mc{realization} = pole(G_q);
    zeros_q_mc{realization} = get_siso_zeros(G_q);

    % Poli e zeri del canale M -> ax
    poles_ax_mc{realization} = pole(G_ax);
    zeros_ax_mc{realization} = get_siso_zeros(G_ax);

    % Salva modello per FRF Monte Carlo
    sys_mc{realization} = sys_ss;

    fprintf('Done\n');

end

%% =========================================================================
%  BLOCCO 2 — CONVERGENZA MONTE CARLO (Task 3.1)

param_names = {'X_u', 'X_q', 'M_u', 'M_q', 'X_\delta', 'M_\delta'};
N_vec = 1:N_realizations;   % asse x: numero di realizzazioni usate
 
% Pre-alloca matrici per media e std cumulative
mean_cumulative = zeros(6, N_realizations);
std_cumulative  = zeros(6, N_realizations);
 
for n = 1:N_realizations
    mean_cumulative(:, n) = mean(theta_mc(:, 1:n), 2);   % media su prime n realiz.
    std_cumulative(:, n)  =  std(theta_mc(:, 1:n), 0, 2); % std su prime n realiz.
end
 
% --- PLOT CONVERGENZA DELLA MEDIA ---
figure('Name', 'MC Convergence: Running Mean', 'Position', [100 100 1100 600]);
for j = 1:6
    subplot(2, 3, j);
    plot(N_vec, mean_cumulative(j, :), 'b-', 'LineWidth', 1.5,'DisplayName', 'Running Mean');
    hold on;
    % Linea orizzontale al valore finale (media su tutte le 50 realizzazioni)
    yline(mean_cumulative(j, end), 'r--', 'LineWidth', 1.2, 'DisplayName', 'Final mean');
    % Banda ±5% attorno al valore finale per valutare la convergenza
    yline(mean_cumulative(j, end) * 1.01, 'k:', 'LineWidth', 1.0, 'HandleVisibility','off');
    yline(mean_cumulative(j, end) * 0.99, 'k:', 'LineWidth', 1.0, 'DisplayName', '±1% band');
    grid on;
    xlabel('Number of realizations N');
    ylabel(sprintf('Mean', j));
    title(sprintf('Running mean: %s', param_names{j}));
    legend('Location', 'best', 'FontSize', 7);
    xlim([1 N_realizations]);
end
sgtitle('Monte Carlo Convergence — Running Mean of Parameters', 'FontSize', 13);
 
% --- PLOT CONVERGENZA DELLA DEVIAZIONE STANDARD ---
figure('Name', 'MC Convergence: Running Std', 'Position', [100 100 1100 600]);
for j = 1:6
    subplot(2, 3, j);
    plot(N_vec, std_cumulative(j, :), 'r-', 'LineWidth', 1.5,'DisplayName', 'Running std');
    hold on;
    yline(std_cumulative(j, end), 'b--', 'LineWidth', 1.2, 'DisplayName', 'Final std');
    yline(std_cumulative(j, end) * 1.05, 'k:', 'LineWidth', 1.0, 'HandleVisibility','off');
    yline(std_cumulative(j, end) * 0.95, 'k:', 'LineWidth', 1.0, 'DisplayName', '±5% band');
    grid on;
    xlabel('Number of realizations N');
    ylabel(sprintf('Std', j));
    title(sprintf('Running std: %s', param_names{j}));
    legend('Location', 'best', 'FontSize', 7);
    xlim([1 N_realizations]);
end
sgtitle('Monte Carlo Convergence — Running Std of Parameters', 'FontSize', 13);
 
% --- STAMPA DIAGNOSTICA: N minimo per convergenza al 5% ---
fprintf('\n========= MC CONVERGENCE ANALYSIS =========\n');
fprintf('Threshold: within 5%% of final value\n\n');
for j = 1:6
    final_mean = mean_cumulative(j, end);
    final_std  = std_cumulative(j, end);
    
    % Trova il primo N oltre il quale la media rimane dentro ±5%
    tol_mean = abs(final_mean) * 0.05;
    conv_mean = find(abs(mean_cumulative(j,:) - final_mean) < tol_mean, 1, 'first');
    
    % Trova il primo N oltre il quale la std rimane dentro ±5%
    tol_std  = abs(final_std) * 0.05;
    conv_std  = find(abs(std_cumulative(j,:)  - final_std)  < tol_std,  1, 'first');
    
    fprintf('  %s:  mean converges at N=%2d,  std converges at N=%2d\n', ...
            param_names{j}, conv_mean, conv_std);
end
fprintf('============================================\n');
%% ========================================================================
%  TASK 3.2 - STATISTICAL ANALYSIS - PARAMETER HISTOGRAMS
% =========================================================================

figure('Name', 'Monte Carlo: Parameter Histograms');

param_names = {'X_u', 'X_q', 'M_u', 'M_q', 'X_\delta', 'M_\delta'};

% Media finale Monte Carlo, cioè media dopo N = 300 realizzazioni
theta_mean_final = mean(theta_mc, 2);

% Errore percentuale tra media Monte Carlo e valore reale
mean_error_abs = theta_mean_final - theta_true;
mean_error_pct = 100 * abs(theta_mean_final - theta_true) ./ abs(theta_true);

% Colore blu scuro per la media
dark_blue = [0.00 0.15 0.55];

for j = 1:6

    subplot(2, 3, j);

    histogram(theta_mc(j, :), 15, ...
    'Normalization', 'pdf', ...
    'FaceColor', [0.2 0.6 0.8], ...
    'EdgeColor', [0.25 0.25 0.25], ...
    'HandleVisibility', 'off');

    hold on;

    % Linea rossa: valore vero
    xline(theta_true(j), 'r--', ...
        'LineWidth', 2, ...
        'DisplayName', 'True value');

    % Linea blu scuro: media Monte Carlo finale, N = 300
    xline(theta_mean_final(j), '--', ...
        'Color', dark_blue, ...
        'LineWidth', 2.2, ...
        'DisplayName', 'MC mean');

    title(sprintf('Distribution of %s', param_names{j}));
    xlabel('Value');
    ylabel('Probability density');
    grid on;

    legend('Location', 'best', 'FontSize', 7);

end

sgtitle('Monte Carlo Parameter Distributions', ...
        'FontSize', 13, 'FontWeight', 'bold');


%% ========================================================================
%  PRINT MEAN ERRORS
% =========================================================================

fprintf('\n========== MONTE CARLO MEAN VS TRUE PARAMETERS ==========\n');
fprintf('Number of realizations: N = %d\n\n', N_realizations);

fprintf('%-12s %-15s %-15s %-15s %-15s\n', ...
    'Parameter', 'True value', 'MC mean', 'Abs error', 'Error [%]');

for j = 1:6

    fprintf('%-12s %-15.6f %-15.6f %-15.6f %-15.3f\n', ...
        param_names{j}, ...
        theta_true(j), ...
        theta_mean_final(j), ...
        mean_error_abs(j), ...
        mean_error_pct(j));

end

fprintf('=========================================================\n');

%% ========================================================================
%  POLE-ZERO MAP: PITCH RATE q
% =========================================================================

figure('Name', 'Poles and Zeros Dispersion - Pitch Rate q', ...
       'Position', [150, 120, 850, 600]);

hold on; grid on; box on;

% Estimated poles
h_poles_q = [];
for j = 1:N_realizations

    p = poles_q_mc{j};

    if ~isempty(p)
        h = plot(real(p), imag(p), ...
            'gx', 'MarkerSize', 8, 'LineWidth', 1.5);

        if isempty(h_poles_q)
            h_poles_q = h;
        end
    end

end

% Estimated zeros
h_zeros_q = [];
for j = 1:N_realizations

    z = zeros_q_mc{j};

    if ~isempty(z)
        h = plot(real(z), imag(z), ...
            'bo', 'MarkerSize', 7, 'LineWidth', 1.5);

        if isempty(h_zeros_q)
            h_zeros_q = h;
        end
    end

end

% True poles
h_true_poles_q = plot(real(poles_true), imag(poles_true), ...
    'rx', 'MarkerSize', 11, 'LineWidth', 2.2);

% True zeros q
h_true_zeros_q = plot(real(zeros_q_true), imag(zeros_q_true), ...
    'ro', 'MarkerSize', 10, 'LineWidth', 2.2);

% Axes
xline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
yline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

xlabel('Real Axis');
ylabel('Imaginary Axis');
title('Poles and Zeros Dispersion in Complex Plane - Pitch Rate q');

legend([h_poles_q, h_zeros_q, h_true_zeros_q, h_true_poles_q], ...
       {'Poles', 'Zeros', 'True Zeros', 'True Poles'}, ...
       'Location', 'best');

axis equal;

%% ========================================================================
%  POLE-ZERO MAP: LONGITUDINAL ACCELERATION ax
% =========================================================================

figure('Name', 'Poles and Zeros Dispersion - Longitudinal Acceleration ax', ...
       'Position', [180, 150, 850, 600]);

hold on; grid on; box on;

% Estimated poles
h_poles_ax = [];
for j = 1:N_realizations

    p = poles_ax_mc{j};

    if ~isempty(p)
        h = plot(real(p), imag(p), ...
            'gx', 'MarkerSize', 8, 'LineWidth', 1.5);

        if isempty(h_poles_ax)
            h_poles_ax = h;
        end
    end

end

% Estimated zeros
h_zeros_ax = [];
for j = 1:N_realizations

    z = zeros_ax_mc{j};

    if ~isempty(z)
        h = plot(real(z), imag(z), ...
            'bo', 'MarkerSize', 7, 'LineWidth', 1.5);

        if isempty(h_zeros_ax)
            h_zeros_ax = h;
        end
    end

end

% True poles
h_true_poles_ax = plot(real(poles_true), imag(poles_true), ...
    'rx', 'MarkerSize', 11, 'LineWidth', 2.2);

% True zeros ax
h_true_zeros_ax = plot(real(zeros_ax_true), imag(zeros_ax_true), ...
    'ro', 'MarkerSize', 10, 'LineWidth', 2.2);

% Axes
xline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
yline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

xlabel('Real Axis');
ylabel('Imaginary Axis');
title('Poles and Zeros Dispersion in Complex Plane - Longitudinal Acceleration a_x');

legend([h_poles_ax, h_zeros_ax, h_true_zeros_ax, h_true_poles_ax], ...
       {'Poles', 'Zeros', 'True Zeros', 'True Poles'}, ...
       'Location', 'best');

axis equal;

%% ========================================================================
%  SORTED POLES/ZEROS FOR HISTOGRAMS
% =========================================================================

% q usually has 2 finite zeros
zeros_q_sorted = NaN(2, N_realizations);

% ax usually has 3 finite zeros
zeros_ax_sorted = NaN(3, N_realizations);

% poles are common: 3 poles
poles_sorted = NaN(3, N_realizations);

for j = 1:N_realizations

    zq = sort_roots_for_hist(zeros_q_mc{j}, 2, 'real_ascending');
    za = sort_roots_for_hist(zeros_ax_mc{j}, 3, 'imag_ascending');
    pp = sort_roots_for_hist(poles_q_mc{j}, 3, 'imag_ascending');

    zeros_q_sorted(:, j)  = zq;
    zeros_ax_sorted(:, j) = za;
    poles_sorted(:, j)    = pp;

end

zeros_q_true_sorted  = sort_roots_for_hist(zeros_q_true, 2, 'real_ascending');
zeros_ax_true_sorted = sort_roots_for_hist(zeros_ax_true, 3, 'imag_ascending');
poles_true_sorted    = sort_roots_for_hist(poles_true, 3, 'imag_ascending');

%% ========================================================================
%  HISTOGRAMS: ZEROS OF PITCH RATE q
% =========================================================================

figure('Name', 'Zeros Histograms - Pitch Rate q', ...
       'Position', [100, 100, 1000, 650]);

sgtitle('Zeros: Pitch Rate q', 'FontSize', 14, 'FontWeight', 'bold');

for k = 1:2

    % Real part
    subplot(2, 2, 2*k - 1);

    data_real = real(zeros_q_sorted(k, :));
    data_real = data_real(~isnan(data_real));

    histogram(data_real, 15, ...
        'FaceColor', [0.2 0.6 0.8], ...
        'EdgeColor', [0.25 0.25 0.25]);

    xline(real(zeros_q_true_sorted(k)), 'r--', 'LineWidth', 2);

    title(sprintf('Distribution Real(Zero %d) - Pitch Rate', k));
    xlabel('Real part');
    ylabel('Count');
    grid on;

    % Imaginary part
    subplot(2, 2, 2*k);

    data_imag = imag(zeros_q_sorted(k, :));
    data_imag = data_imag(~isnan(data_imag));

    histogram(data_imag, 15, ...
        'FaceColor', [0.2 0.6 0.8], ...
        'EdgeColor', [0.25 0.25 0.25]);

    xline(imag(zeros_q_true_sorted(k)), 'r--', 'LineWidth', 2);

    title(sprintf('Distribution Imag(Zero %d) - Pitch Rate', k));
    xlabel('Imaginary part');
    ylabel('Count');
    grid on;

end

%% ========================================================================
%  HISTOGRAMS: ZEROS OF LONGITUDINAL ACCELERATION ax
% =========================================================================

figure('Name', 'Zeros Histograms - Longitudinal Acceleration ax', ...
       'Position', [100, 100, 1050, 750]);

sgtitle('Zeros: Longitudinal Acceleration a_x', ...
        'FontSize', 14, 'FontWeight', 'bold');

for k = 1:3

    % Real part
    subplot(3, 2, 2*k - 1);

    data_real = real(zeros_ax_sorted(k, :));
    data_real = data_real(~isnan(data_real));

    histogram(data_real, 15, ...
        'FaceColor', [0.2 0.6 0.8], ...
        'EdgeColor', [0.25 0.25 0.25]);

    xline(real(zeros_ax_true_sorted(k)), 'r--', 'LineWidth', 2);

    title(sprintf('Distribution Real(Zero %d) - Longitudinal Acceleration', k));
    xlabel('Real part');
    ylabel('Count');
    grid on;

    % Imaginary part
    subplot(3, 2, 2*k);

    data_imag = imag(zeros_ax_sorted(k, :));
    data_imag = data_imag(~isnan(data_imag));

    histogram(data_imag, 15, ...
        'FaceColor', [0.2 0.6 0.8], ...
        'EdgeColor', [0.25 0.25 0.25]);

    xline(imag(zeros_ax_true_sorted(k)), 'r--', 'LineWidth', 2);

    title(sprintf('Distribution Imag(Zero %d) - Longitudinal Acceleration', k));
    xlabel('Imaginary part');
    ylabel('Count');
    grid on;

end

%% ========================================================================
%  HISTOGRAMS: POLES
% =========================================================================

figure('Name', 'Poles Histograms', ...
       'Position', [100, 100, 1200, 650]);

sgtitle('Poles: Pitch Rate and Longitudinal Acceleration', ...
        'FontSize', 14, 'FontWeight', 'bold');

for k = 1:3

    % Real part
    subplot(2, 3, k);

    data_real = real(poles_sorted(k, :));
    data_real = data_real(~isnan(data_real));

    histogram(data_real, 15, ...
        'FaceColor', [0.2 0.6 0.8], ...
        'EdgeColor', [0.25 0.25 0.25]);

    xline(real(poles_true_sorted(k)), 'r--', 'LineWidth', 2);

    title(sprintf('Poles Real Part - Pole %d', k));
    xlabel('Real part');
    ylabel('Count');
    grid on;

    % Imaginary part
    subplot(2, 3, k + 3);

    data_imag = imag(poles_sorted(k, :));
    data_imag = data_imag(~isnan(data_imag));

    histogram(data_imag, 15, ...
        'FaceColor', [0.2 0.6 0.8], ...
        'EdgeColor', [0.25 0.25 0.25]);

    xline(imag(poles_true_sorted(k)), 'r--', 'LineWidth', 2);

    title(sprintf('Poles Imaginary Part - Pole %d', k));
    xlabel('Imaginary part');
    ylabel('Count');
    grid on;

end

%% ========================================================================
%  PRINT TRUE POLES AND ZEROS
% =========================================================================

fprintf('\n========== TRUE POLES AND ZEROS ==========\n');

fprintf('\nTrue poles:\n');
disp(poles_true_sorted);

fprintf('\nTrue zeros q:\n');
disp(zeros_q_true_sorted);

fprintf('\nTrue zeros ax:\n');
disp(zeros_ax_true_sorted);

%% ========================================================================
%  MONTE CARLO FRF CLOUD - STYLE LIKE REFERENCE
% =========================================================================

freq_plot_Hz = logspace(-1, 1, 400);
w_plot = 2*pi*freq_plot_Hz;

mag_q_all    = zeros(N_realizations, length(w_plot));
phase_q_all  = zeros(N_realizations, length(w_plot));

mag_ax_all   = zeros(N_realizations, length(w_plot));
phase_ax_all = zeros(N_realizations, length(w_plot));

for j = 1:N_realizations

    Gq_j  = sys_mc{j}(1,1);
    Gax_j = sys_mc{j}(2,1);

    Hq  = squeeze(freqresp(Gq_j,  w_plot));
    Hax = squeeze(freqresp(Gax_j, w_plot));

    mag_q_all(j, :)    = 20*log10(abs(Hq));
    phase_q_all(j, :)  = unwrap(angle(Hq))*180/pi;

    mag_ax_all(j, :)   = 20*log10(abs(Hax));
    phase_ax_all(j, :) = unwrap(angle(Hax))*180/pi;

end

% True FRF
Hq_true  = squeeze(freqresp(Gq_true,  w_plot));
Hax_true = squeeze(freqresp(Gax_true, w_plot));

mag_q_true    = 20*log10(abs(Hq_true));
phase_q_true  = unwrap(angle(Hq_true))*180/pi;

mag_ax_true   = 20*log10(abs(Hax_true));
phase_ax_true = unwrap(angle(Hax_true))*180/pi;

% Mean and std
mag_q_mean    = mean(mag_q_all, 1);
phase_q_mean  = mean(phase_q_all, 1);

mag_ax_mean   = mean(mag_ax_all, 1);
phase_ax_mean = mean(phase_ax_all, 1);

mag_q_std     = std(mag_q_all, 0, 1);
phase_q_std   = std(phase_q_all, 0, 1);

mag_ax_std    = std(mag_ax_all, 0, 1);
phase_ax_std  = std(phase_ax_all, 0, 1);

figure('Name', 'Monte Carlo FRF Cloud', ...
       'Position', [100, 80, 1300, 750]);

% ------------------------------------------------------------------------
% Pitch rate magnitude
% ------------------------------------------------------------------------
subplot(2,2,1);
hold on; grid on; box on;

for j = 1:N_realizations
    semilogx(freq_plot_Hz, mag_q_all(j,:), ...
        'Color', [0.75 0.75 0.75], 'LineWidth', 0.6, ...
        'HandleVisibility', 'off');
end

h_mean = semilogx(freq_plot_Hz, mag_q_mean, ...
    'r-', 'LineWidth', 2);

h_true = semilogx(freq_plot_Hz, mag_q_true, ...
    'b--', 'LineWidth', 2);

h_3s = semilogx(freq_plot_Hz, mag_q_mean + 3*mag_q_std, ...
    'k:', 'LineWidth', 1.5);

semilogx(freq_plot_Hz, mag_q_mean - 3*mag_q_std, ...
    'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');

title('Dispersion of Pitch Rate FRF');
ylabel('Magnitude [dB]');
legend([h_mean, h_true, h_3s], ...
       {'Mean', 'True Model', '\pm3\sigma Bounds'}, ...
       'Location', 'best');

% ------------------------------------------------------------------------
% Longitudinal acceleration magnitude
% ------------------------------------------------------------------------
subplot(2,2,2);
hold on; grid on; box on;

for j = 1:N_realizations
    semilogx(freq_plot_Hz, mag_ax_all(j,:), ...
        'Color', [0.75 0.75 0.75], 'LineWidth', 0.6, ...
        'HandleVisibility', 'off');
end

h_mean = semilogx(freq_plot_Hz, mag_ax_mean, ...
    'r-', 'LineWidth', 2);

h_true = semilogx(freq_plot_Hz, mag_ax_true, ...
    'b--', 'LineWidth', 2);

h_3s = semilogx(freq_plot_Hz, mag_ax_mean + 3*mag_ax_std, ...
    'k:', 'LineWidth', 1.5);

semilogx(freq_plot_Hz, mag_ax_mean - 3*mag_ax_std, ...
    'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');

title('Dispersion of Longitudinal Acceleration FRF');
ylabel('Magnitude [dB]');
legend([h_mean, h_true, h_3s], ...
       {'Mean', 'True Model', '\pm3\sigma Bounds'}, ...
       'Location', 'best');

% ------------------------------------------------------------------------
% Pitch rate phase
% ------------------------------------------------------------------------
subplot(2,2,3);
hold on; grid on; box on;

for j = 1:N_realizations
    semilogx(freq_plot_Hz, phase_q_all(j,:), ...
        'Color', [0.75 0.75 0.75], 'LineWidth', 0.6, ...
        'HandleVisibility', 'off');
end

semilogx(freq_plot_Hz, phase_q_mean, ...
    'r-', 'LineWidth', 2);

semilogx(freq_plot_Hz, phase_q_true, ...
    'b--', 'LineWidth', 2);

semilogx(freq_plot_Hz, phase_q_mean + 3*phase_q_std, ...
    'k:', 'LineWidth', 1.5);

semilogx(freq_plot_Hz, phase_q_mean - 3*phase_q_std, ...
    'k:', 'LineWidth', 1.5);

xlabel('Frequency [Hz]');
ylabel('Phase [deg]');

% ------------------------------------------------------------------------
% Longitudinal acceleration phase
% ------------------------------------------------------------------------
subplot(2,2,4);
hold on; grid on; box on;

for j = 1:N_realizations
    semilogx(freq_plot_Hz, phase_ax_all(j,:), ...
        'Color', [0.75 0.75 0.75], 'LineWidth', 0.6, ...
        'HandleVisibility', 'off');
end

semilogx(freq_plot_Hz, phase_ax_mean, ...
    'r-', 'LineWidth', 2);

semilogx(freq_plot_Hz, phase_ax_true, ...
    'b--', 'LineWidth', 2);

semilogx(freq_plot_Hz, phase_ax_mean + 3*phase_ax_std, ...
    'k:', 'LineWidth', 1.5);

semilogx(freq_plot_Hz, phase_ax_mean - 3*phase_ax_std, ...
    'k:', 'LineWidth', 1.5);

xlabel('Frequency [Hz]');
ylabel('Phase [deg]');

%% ========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function [A, B, C, D] = System_matrix(theta, T)

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

    C = [0,   1, 0;
         X_u, X_q, 0];

    D = [0;
         X_d];

end

function z = get_siso_zeros(G)

    % Converte la singola FRF in transfer function
    G_tf = tf(G);

    % Estrae numeratore
    [num, ~] = tfdata(G_tf, 'v');

    % Tolleranza numerica
    tol = 1e-9;

    if isempty(num)
        z = [];
        return
    end

    max_num = max(abs(num));

    if max_num == 0 || isnan(max_num)
        z = [];
        return
    end

    % Elimina coefficienti numericamente nulli
    num(abs(num) < tol * max_num) = 0;

    % Rimuove zeri iniziali nel polinomio
    first_nonzero = find(abs(num) > 0, 1, 'first');

    if isempty(first_nonzero)
        z = [];
        return
    end

    num = num(first_nonzero:end);

    % Se il numeratore è costante, non ci sono zeri finiti
    if length(num) <= 1
        z = [];
        return
    end

    % Zeri come radici del numeratore
    z = roots(num);

    % Rimuove eventuali NaN/Inf
    z = z(isfinite(real(z)) & isfinite(imag(z)));

end

function z_sorted = sort_roots_for_hist(z, n_expected, mode)

    z_sorted = NaN(n_expected, 1);

    if isempty(z)
        return
    end

    z = z(:);
    z = z(isfinite(real(z)) & isfinite(imag(z)));

    if isempty(z)
        return
    end

    switch mode

        case 'real_ascending'
            [~, idx] = sort(real(z), 'ascend');
            z = z(idx);

        case 'imag_ascending'
            [~, idx] = sortrows([imag(z), real(z)], [1 2]);
            z = z(idx);

        otherwise
            [~, idx] = sort(real(z), 'ascend');
            z = z(idx);

    end

    n_copy = min(length(z), n_expected);
    z_sorted(1:n_copy) = z(1:n_copy);

end