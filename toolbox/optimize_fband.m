function result = optimize_fband(u, q, ax, Ts, gamma_sq_q, gamma_sq_ax, f_coh, init_sys, opt)
    % ============================================================
    % DISABILITAZIONE WARNING
    % ============================================================
    orig_warning_state = warning('query', 'all');
    warning('off', 'all');

    % ============================================================
    % Parametri di configurazione
    % ============================================================
    coh_thr    = 0.5;   % soglia minima di coerenza per definire la banda base
    step_Hz    = 0.2;   % passo della griglia di ricerca [Hz]
    min_points = 8;     % numero minimo di punti frequenziali per banda valida
    search_range = 3.0; % [Hz] quanto espandere la ricerca attorno alla banda base

    % ============================================================
    % 1) Banda iniziale da coerenza
    %    Trova il blocco contiguo più lungo in cui ENTRAMBE le
    %    coerenze superano coh_thr — usato come punto di partenza
    %    per la grid search.
    % ============================================================
    valid = (gamma_sq_q > coh_thr) & (gamma_sq_ax > coh_thr);
    if ~any(valid)
        warning(orig_warning_state);
        error('Non esistono frequenze in cui entrambe le coherence sono > %.2f', coh_thr);
    end
    d = diff([false; valid; false]);
    start_idx    = find(d ==  1);
    end_idx      = find(d == -1) - 1;
    block_lengths = end_idx - start_idx + 1;
    [~, k_best]  = max(block_lengths);
    f_min_0 = f_coh(start_idx(k_best));
    f_max_0 = f_coh(end_idx(k_best));
    fprintf('\nBanda base (coherence > %.2f): [%.2f, %.2f] Hz\n', coh_thr, f_min_0, f_max_0);

    % ============================================================
    % 2) Pre-calcolo FFT su tutti i segnali (una volta sola)
    % ============================================================
    N    = length(u);
    fs   = 1 / Ts;
    f_nyq = fs / 2;
    U    = fft(u);
    Y    = fft([q, ax]);
    Npos = floor(N/2) + 1;
    U_pos   = U(1:Npos);
    Y_pos   = Y(1:Npos, :);
    freq_Hz = (0:Npos-1)' * fs / N;
    W_pos   = 2*pi * freq_Hz;

    % ============================================================
    % 3) Costruzione griglia e pre-filtraggio delle combinazioni
    %    valide (almeno min_points punti FFT nella banda)
    % ============================================================
    f_min_vec = max(0.2, f_min_0 - search_range) : step_Hz : f_min_0;
    f_max_vec = f_max_0 : step_Hz : min(f_nyq, f_max_0 + search_range);
    [F_MIN, F_MAX] = meshgrid(f_min_vec, f_max_vec);
    f_min_raw = F_MIN(:);
    f_max_raw = F_MAX(:);
    num_raw   = length(f_min_raw);

    valid_mask = (f_min_raw < f_max_raw);
    for k = 1:num_raw
        if valid_mask(k)
            n_pts = sum((freq_Hz >= f_min_raw(k)) & (freq_Hz <= f_max_raw(k)));
            if n_pts < min_points
                valid_mask(k) = false;
            end
        end
    end

    valid_indices    = find(valid_mask);
    f_min_par        = f_min_raw(valid_indices);
    f_max_par        = f_max_raw(valid_indices);
    num_combinations = length(f_min_par);

    if num_combinations == 0
        warning(orig_warning_state);
        error('Nessuna combinazione (f_min, f_max) soddisfa i criteri minimi.');
    end

    % ============================================================
    % Array pre-allocati per il ciclo parfor
    % PASSO 1: raccolta valori GREZZI (weighted_fit e mean_rel_std)
    %          La normalizzazione z-score avviene DOPO il parfor,
    %          quando tutti i valori sono disponibili.
    % ============================================================
    wf_par      = -Inf(num_combinations, 1);  % weighted fit grezzo
    rel_std_par = -Inf(num_combinations, 1);  % deviazione std relativa media grezza
    fit_q_par   = -Inf(num_combinations, 1);
    fit_ax_par  = -Inf(num_combinations, 1);

    fprintf('Inizio grid search: %d combinazioni saltate, %d in elaborazione parfor...\n', ...
        num_raw - num_combinations, num_combinations);

    % ============================================================
    % DataQueue per mostrare avanzamento parfor
    % ============================================================
    dq    = parallel.pool.DataQueue;
    count = 0;
    afterEach(dq, @updateProgress);
    function updateProgress(~)
        count      = count + 1;
        step_print = max(1, round(num_combinations / 20));
        if mod(count, step_print) == 0 || count == num_combinations
            fprintf('Progresso: %d / %d (%.1f%%)\n', count, num_combinations, ...
                    100 * count / num_combinations);
        end
    end

    % ============================================================
    % 4) Ciclo parallelo — SOLO raccolta valori grezzi
    % ============================================================
    opt_parfor         = opt;
    opt_parfor.Display = 'off';

    parfor i = 1:num_combinations
        warning('off', 'all');

        f_min_curr = f_min_par(i);
        f_max_curr = f_max_par(i);

        % Estrazione dati in banda
        idx = (freq_Hz >= f_min_curr) & (freq_Hz <= f_max_curr);
        U_f = U_pos(idx);
        Y_f = Y_pos(idx, :);
        W_f = W_pos(idx);

        % iddata frequenziale (Ts=0: non usato con 'Frequency' esplicito)
        data = iddata(Y_f, U_f, 0, 'Frequency', W_f);

        try
            sys_est = greyest(data, init_sys, opt_parfor);

            fit_q_i  = sys_est.Report.Fit.FitPercent(1);
            fit_ax_i = sys_est.Report.Fit.FitPercent(2);

            % --- Coerenza media nella banda ---
            idx_coh     = (f_coh >= f_min_curr) & (f_coh <= f_max_curr);
            if any(idx_coh)
                mean_coh_q  = mean(gamma_sq_q(idx_coh));
                mean_coh_ax = mean(gamma_sq_ax(idx_coh));
            else
                mean_coh_q  = eps;
                mean_coh_ax = eps;
            end

            % --- Fit pesato per coerenza ---
            sum_coh = mean_coh_q + mean_coh_ax;
            if sum_coh > 0
                wf_i = (fit_q_i * mean_coh_q + fit_ax_i * mean_coh_ax) / sum_coh;
            else
                wf_i = eps;
            end

            % --- Deviazione standard relativa media dei parametri ---
            theta_i   = getpvec(sys_est, 'free');
            cov_mat_i = getcov(sys_est);

            if isempty(cov_mat_i) || any(isnan(cov_mat_i(:))) || any(isinf(cov_mat_i(:)))
                rel_std_i = Inf;
            else
                std_i        = sqrt(diag(cov_mat_i));
                theta_safe_i = theta_i;
                theta_safe_i(theta_safe_i == 0) = eps;
                rel_std_i = mean(100 * std_i ./ abs(theta_safe_i));
            end

            % Salva solo se i valori sono finiti e validi
            if isfinite(wf_i) && isfinite(rel_std_i) && rel_std_i > 0
                wf_par(i)      = wf_i;
                rel_std_par(i) = rel_std_i;
                fit_q_par(i)   = fit_q_i;
                fit_ax_par(i)  = fit_ax_i;
            end

        catch
            % Lascia i valori a -Inf in caso di fallimento di greyest
        end

        send(dq, i);
    end

    % ============================================================
    % 5) Normalizzazione z-score POST parfor
    %    Solo ora abbiamo tutti i valori e possiamo calcolare
    %    media e std sull'intera griglia.
    % ============================================================
    valid_results = (wf_par > -Inf) & (rel_std_par > -Inf);

    if sum(valid_results) < 2
        warning(orig_warning_state);
        error(['Meno di 2 combinazioni valide (%d trovate). ' ...
               'Verificare theta_0, banda di coerenza o min_points.'], ...
               sum(valid_results));
    end

    fprintf('\n%d / %d combinazioni hanno prodotto una stima valida.\n', ...
            sum(valid_results), num_combinations);

    % Statistiche sui soli valori validi
    mu_wf  = mean(wf_par(valid_results));
    s_wf   =  std(wf_par(valid_results));
    mu_std = mean(rel_std_par(valid_results));
    s_std  =  std(rel_std_par(valid_results));

    % Gestione caso degenere (tutti i valori identici -> std = 0)
    if s_wf  < eps, s_wf  = 1; end
    if s_std < eps, s_std = 1; end

    % Calcolo score normalizzato per ogni combinazione valida
    score_par = -Inf(num_combinations, 1);
    for i = 1:num_combinations
        if valid_results(i)
            wf_norm  = (wf_par(i)      - mu_wf)  / s_wf;
            std_norm = (rel_std_par(i) - mu_std) / s_std;
            % Massimizza fit normalizzato, minimizza incertezza normalizzata
            score_par(i) = wf_norm - std_norm;
        end
    end

    % ============================================================
    % 6) Ricostruzione griglia originale e individuazione massimo
    % ============================================================
    score_flat = -Inf(num_raw, 1);
    score_flat(valid_indices) = score_par;

    fit_q_flat = -Inf(num_raw, 1);
    fit_q_flat(valid_indices) = fit_q_par;

    fit_ax_flat = -Inf(num_raw, 1);
    fit_ax_flat(valid_indices) = fit_ax_par;

    [best_score, best_idx] = max(score_flat);
    if best_score == -Inf
        warning(orig_warning_state);
        error('Nessuna combinazione ha prodotto uno score valido.');
    end

    best_f_min = f_min_raw(best_idx);
    best_f_max = f_max_raw(best_idx);

    % Output
    result.f_min    = best_f_min;
    result.f_max    = best_f_max;
    result.fit_q    = fit_q_flat(best_idx);
    result.fit_ax   = fit_ax_flat(best_idx);
    result.fit_sum  = fit_q_flat(best_idx) + fit_ax_flat(best_idx);
    result.score    = best_score;

    fprintf('\n========= BANDA OTTIMA TROVATA =========\n');
    fprintf('f_min   = %.2f Hz\n', result.f_min);
    fprintf('f_max   = %.2f Hz\n', result.f_max);
    fprintf('fit_q   = %.2f %%\n', result.fit_q);
    fprintf('fit_ax  = %.2f %%\n', result.fit_ax);
    fprintf('Score   = %.4f\n',    result.score);
    fprintf('=========================================\n');

    % ============================================================
    % 7) Plot 3D score normalizzato
    % ============================================================
    SCORE_MATRIX = reshape(score_flat, size(F_MIN));
    SCORE_MATRIX(SCORE_MATRIX == -Inf) = NaN;

    figure('Name', 'Ottimizzazione f_min e f_max', 'NumberTitle', 'off');
    surf(F_MIN, F_MAX, SCORE_MATRIX, 'EdgeColor', 'interp', 'FaceColor', 'interp');
    hold on;
    plot3(best_f_min, best_f_max, best_score, 'r*', 'MarkerSize', 15, 'LineWidth', 2);
    title('Grid Search: Score z-normalizzato (fit pesato - rel. std) vs Frequenze');
    xlabel('f_{min} (Hz)');
    ylabel('f_{max} (Hz)');
    zlabel('Score (adim.)');
    colorbar;
    view(-45, 35);
    grid on;

    warning(orig_warning_state);
end