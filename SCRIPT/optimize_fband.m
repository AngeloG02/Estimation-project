
function result = optimize_fband(u, q, ax, Ts, gamma_sq_q, gamma_sq_ax, f_coh, init_sys, opt)

    % Parametri
    coh_thr = 0.5;
    step_Hz = 0.1;
    min_points = 8;   % minimo numero di punti in frequenza per evitare stime assurde

    % ============================================================
    % 1) Banda iniziale: frequenze dove entrambe le coherence > 0.5
    % ============================================================
    valid = (gamma_sq_q > coh_thr) & (gamma_sq_ax > coh_thr);

if ~any(valid)
    error('Non esistono frequenze in cui entrambe le coherence sono > %.2f', coh_thr);
end

% Trova i blocchi continui di frequenze valide
d = diff([false; valid; false]);
start_idx = find(d == 1);
end_idx   = find(d == -1) - 1;

% Scegli il blocco continuo più lungo
block_lengths = end_idx - start_idx + 1;
[~, k_best] = max(block_lengths);

f_min_0 = f_coh(start_idx(k_best));
f_max_0 = f_coh(end_idx(k_best));

    fprintf('\nBanda iniziale da coherence > %.2f: [%.2f, %.2f] Hz\n', coh_thr, f_min_0, f_max_0);

    % ============================================================
    % 2) Valutazione iniziale
    % ============================================================
    [sys_best, fit_q_best, fit_ax_best, fit_sum_best, ok] = ...
        estimate_on_band(u, q, ax, Ts, f_min_0, f_max_0, init_sys, opt, min_points);

    if ~ok
        error('La stima iniziale sulla banda [%.2f, %.2f] Hz è fallita.', f_min_0, f_max_0);
    end

    f_min_best = f_min_0;
    f_max_best = f_max_0;

    fprintf('Banda iniziale -> fit_q = %.2f, fit_ax = %.2f, somma = %.2f\n', ...
        fit_q_best, fit_ax_best, fit_sum_best);

    % ============================================================
    % 3) Ottimizzazione di f_min: scendo a step di 0.2 Hz
    % ============================================================
    improved = true;
    while improved
        improved = false;
        f_min_candidate = f_min_best - step_Hz;

        if f_min_candidate <= 0
            break;
        end

        [sys_cand, fit_q_cand, fit_ax_cand, fit_sum_cand, ok] = ...
            estimate_on_band(u, q, ax, Ts, f_min_candidate, f_max_best, init_sys, opt, min_points);

        if ok && (fit_sum_cand > fit_sum_best)
            fprintf('Migliora abbassando f_min: %.2f -> %.2f | somma fit %.2f -> %.2f\n', ...
                f_min_best, f_min_candidate, fit_sum_best, fit_sum_cand);

            f_min_best = f_min_candidate;
            sys_best = sys_cand;
            fit_q_best = fit_q_cand;
            fit_ax_best = fit_ax_cand;
            fit_sum_best = fit_sum_cand;
            improved = true;
        end
    end

    % ============================================================
    % 4) Ottimizzazione di f_max: salgo a step di 0.2 Hz
    % ============================================================
    N = length(u);
    fs = 1/Ts;
    f_nyq = fs/2;

    improved = true;
    while improved
        improved = false;
        f_max_candidate = f_max_best + step_Hz;

        if f_max_candidate >= f_nyq
            break;
        end

        [sys_cand, fit_q_cand, fit_ax_cand, fit_sum_cand, ok] = ...
            estimate_on_band(u, q, ax, Ts, f_min_best, f_max_candidate, init_sys, opt, min_points);

        if ok && (fit_sum_cand > fit_sum_best)
            fprintf('Migliora alzando f_max: %.2f -> %.2f | somma fit %.2f -> %.2f\n', ...
                f_max_best, f_max_candidate, fit_sum_best, fit_sum_cand);

            f_max_best = f_max_candidate;
            sys_best = sys_cand;
            fit_q_best = fit_q_cand;
            fit_ax_best = fit_ax_cand;
            fit_sum_best = fit_sum_cand;
            improved = true;
        end
    end

    % ============================================================
    % 5) Output finale
    % ============================================================
    result.f_min = f_min_best;
    result.f_max = f_max_best;
    result.sys = sys_best;
    result.fit_q = fit_q_best;
    result.fit_ax = fit_ax_best;
    result.fit_sum = fit_sum_best;

    fprintf('\n========= BANDA OTTIMA TROVATA =========\n');
    fprintf('f_min = %.2f Hz\n', result.f_min);
    fprintf('f_max = %.2f Hz\n', result.f_max);
    fprintf('fit_q = %.2f %%\n', result.fit_q);
    fprintf('fit_ax = %.2f %%\n', result.fit_ax);
    fprintf('fit_sum = %.2f %%\n', result.fit_sum);
end


function [sys, fit_q, fit_ax, fit_sum, ok] = estimate_on_band(u, q, ax, Ts, f_min, f_max, init_sys, opt, min_points)

    ok = false;
    sys = [];
    fit_q = -Inf;
    fit_ax = -Inf;
    fit_sum = -Inf;

    try
        % ============================================================
        % 1) FFT solo su frequenze positive
        % ============================================================
        N = length(u);
        fs = 1/Ts;

        U = fft(u);
        Y = fft([q, ax]);

        Npos = floor(N/2) + 1;
        U = U(1:Npos);
        Y = Y(1:Npos, :);

        freq_Hz = (0:Npos-1)' * fs / N;
        W = 2*pi*freq_Hz;

        % ============================================================
        % 2) Selezione banda
        % ============================================================
        idx = (freq_Hz >= f_min) & (freq_Hz <= f_max);

        if sum(idx) < min_points
            fprintf('Banda [%.2f, %.2f] scartata: troppo pochi punti (%d)\n', ...
                f_min, f_max, sum(idx));
            return;
        end

        U_f = U(idx);
        Y_f = Y(idx, :);
        W_f = W(idx);

        % ============================================================
        % 3) Crea iddata frequenziale
        % ============================================================
        data = iddata(Y_f, U_f, 0, 'Frequency', W_f);
        data.InputName = 'Pitch Moment';
        data.OutputName = {'Pitch rate', 'Longitudinal acceleration'};
        data.OutputUnit = {'rad/s', 'm/s^2'};
        data.TimeUnit = 's';

        % ============================================================
        % 4) Stima grey-box
        % ============================================================
        sys = greyest(data, init_sys, opt);

        % ============================================================
        % 5) Calcolo risposta del modello sulle stesse frequenze
        % ============================================================
        % G = freqresp(sys, W_f);   % dimensione: ny x nu x Nf
        % G = squeeze(G);           % diventa ny x Nf se nu = 1
        % 
        % if size(G,1) ~= 2
        %     error('La dimensione della FRF del modello non è coerente con 2 uscite.');
        % end

        % Yhat_q  = (G(1,:).') .* U_f;
        % Yhat_ax = (G(2,:).') .* U_f;
        % 
        % Ytrue_q  = Y_f(:,1);
        % Ytrue_ax = Y_f(:,2);
        % 
        % % ============================================================
        % % 6) Fit percentuale
        % % ============================================================
        % fit_q = 100 * (1 - norm(Ytrue_q - Yhat_q) / norm(Ytrue_q - mean(Ytrue_q)));
        % fit_ax = 100 * (1 - norm(Ytrue_ax - Yhat_ax) / norm(Ytrue_ax - mean(Ytrue_ax)));
        
        fit_sum = sum(sys.Report.Fit.FitPercent);
        fit_ax  = sys.Report.Fit.FitPercent(2);
        fit_q = sys.Report.Fit.FitPercent(1);
        ok = true;

    catch ME
        fprintf('Stima fallita su banda [%.2f, %.2f] Hz -> %s\n', f_min, f_max, ME.message);
    end
end