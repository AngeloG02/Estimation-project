function result = optimize_fband_grid(u, q, ax, Ts, gamma_sq_q, gamma_sq_ax, f_coh, init_sys, opt)

    % ============================================================
    % DISABILITAZIONE WARNING
    % ============================================================
    % Salviamo lo stato corrente per ripristinarlo alla fine
    orig_warning_state = warning('query', 'all');
    warning('off', 'all');

    % ============================================================
    % Parametri di configurazione
    % ============================================================
    coh_thr = 0.5;
    step_Hz = 0.1;
    min_points = 8;     % Minimo numero di punti in frequenza validi
    search_range = 2.0; % +/- 2 Hz dalla frequenza iniziale

    % ============================================================
    % 1) Banda iniziale: frequenze dove entrambe le coherence > 0.5
    % ============================================================
    valid = (gamma_sq_q > coh_thr) & (gamma_sq_ax > coh_thr);

    if ~any(valid)
        warning(orig_warning_state); % Ripristina prima dell'errore
        error('Non esistono frequenze in cui entrambe le coherence sono > %.2f', coh_thr);
    end

    d = diff([false; valid; false]);
    start_idx = find(d == 1);
    end_idx   = find(d == -1) - 1;

    block_lengths = end_idx - start_idx + 1;
    [~, k_best] = max(block_lengths);

    f_min_0 = f_coh(start_idx(k_best));
    f_max_0 = f_coh(end_idx(k_best));

    fprintf('\nBanda base (coherence > %.2f): [%.2f, %.2f] Hz\n', coh_thr, f_min_0, f_max_0);

    % ============================================================
    % 2) Pre-calcolo FFT su tutti i segnali (OTTIMIZZAZIONE)
    % ============================================================
    N = length(u);
    fs = 1/Ts;
    f_nyq = fs/2;

    U = fft(u);
    Y = fft([q, ax]);

    Npos = floor(N/2) + 1;
    U_pos = U(1:Npos);
    Y_pos = Y(1:Npos, :);

    freq_Hz = (0:Npos-1)' * fs / N;
    W_pos = 2*pi*freq_Hz;

    % ============================================================
    % 3) Costruzione della griglia di ricerca
    % ============================================================
    f_min_vec = max(0.01, f_min_0 - search_range) : step_Hz : f_min_0;
    f_max_vec = f_max_0 : step_Hz : min(f_nyq, f_max_0 + search_range);

    [F_MIN, F_MAX] = meshgrid(f_min_vec, f_max_vec);
    f_min_flat = F_MIN(:);
    f_max_flat = F_MAX(:);
    num_combinations = length(f_min_flat);

    fit_sum_flat = -Inf(num_combinations, 1);
    fit_q_flat   = -Inf(num_combinations, 1);
    fit_ax_flat  = -Inf(num_combinations, 1);

    fprintf('Inizio calcolo a griglia su %d combinazioni con parfor...\n', num_combinations);

    % ============================================================
    % SETUP: DataQueue per mostrare l'avanzamento del parfor
    % ============================================================
    dq = parallel.pool.DataQueue;
    count = 0;
    afterEach(dq, @updateProgress);

    function updateProgress(~)
        count = count + 1;
        step_print = max(1, round(num_combinations / 20)); 
        if mod(count, step_print) == 0 || count == num_combinations
            fprintf('Progresso: %d / %d (%.1f%%)\n', count, num_combinations, (count/num_combinations)*100);
        end
    end

    % ============================================================
    % 4) Ciclo Parallelo (Grid Search)
    % ============================================================
    opt_parfor = opt;
    opt_parfor.Display = 'off';

    parfor i = 1:num_combinations
        % Assicuriamoci che anche i worker paralleli ignorino i warning
        warning('off', 'all'); 
        
        f_min_curr = f_min_flat(i);
        f_max_curr = f_max_flat(i);

        if f_min_curr >= f_max_curr
            send(dq, i); 
            continue;
        end

        idx = (freq_Hz >= f_min_curr) & (freq_Hz <= f_max_curr);

        if sum(idx) >= min_points
            U_f = U_pos(idx);
            Y_f = Y_pos(idx, :);
            W_f = W_pos(idx);

            data = iddata(Y_f, U_f, 0, 'Frequency', W_f);
            data.InputName = 'Pitch Moment';
            data.OutputName = {'Pitch rate', 'Longitudinal acceleration'};
            data.OutputUnit = {'rad/s', 'm/s^2'};
            data.TimeUnit = 's';

            try
                sys_est = greyest(data, init_sys, opt_parfor);
                
                fit_sum_flat(i) = sum(sys_est.Report.Fit.FitPercent);
                fit_q_flat(i)   = sys_est.Report.Fit.FitPercent(1);
                fit_ax_flat(i)  = sys_est.Report.Fit.FitPercent(2);
            catch
                % Fallito, lascia a -Inf
            end
        end
        
        send(dq, i); 
    end

    % ============================================================
    % 5) Individuazione del massimo globale
    % ============================================================
    [best_fit_sum, best_idx] = max(fit_sum_flat);

    if best_fit_sum == -Inf
        warning(orig_warning_state); % Ripristina
        error('Nessuna combinazione ha prodotto una stima valida.');
    end

    best_f_min = f_min_flat(best_idx);
    best_f_max = f_max_flat(best_idx);

    % ============================================================
    % 6) Calcolo del modello finale sulla combinazione ottima
    % ============================================================
    idx_best = (freq_Hz >= best_f_min) & (freq_Hz <= best_f_max);
    data_best = iddata(Y_pos(idx_best, :), U_pos(idx_best), 0, 'Frequency', W_pos(idx_best));
    data_best.InputName = 'Pitch Moment';
    data_best.OutputName = {'Pitch rate', 'Longitudinal acceleration'};
    data_best.OutputUnit = {'rad/s', 'm/s^2'};
    data_best.TimeUnit = 's';
    
    sys_best = greyest(data_best, init_sys, opt);

    result.f_min = best_f_min;
    result.f_max = best_f_max;
    result.sys = sys_best;
    result.fit_q = fit_q_flat(best_idx);
    result.fit_ax = fit_ax_flat(best_idx);
    result.fit_sum = best_fit_sum;

    fprintf('\n========= BANDA OTTIMA TROVATA =========\n');
    fprintf('f_min   = %.2f Hz\n', result.f_min);
    fprintf('f_max   = %.2f Hz\n', result.f_max);
    fprintf('fit_q   = %.2f %%\n', result.fit_q);
    fprintf('fit_ax  = %.2f %%\n', result.fit_ax);
    fprintf('fit_sum = %.2f %%\n', result.fit_sum);

    % ============================================================
    % 7) Plot 3D dei risultati
    % ============================================================
    FIT_SUM_MATRIX = reshape(fit_sum_flat, size(F_MIN));
    FIT_SUM_MATRIX(FIT_SUM_MATRIX == -Inf) = NaN; 

    figure('Name', 'Ottimizzazione f_min e f_max', 'NumberTitle', 'off');
    surf(F_MIN, F_MAX, FIT_SUM_MATRIX, 'EdgeColor', 'interp', 'FaceColor', 'interp');
    hold on;
    plot3(best_f_min, best_f_max, best_fit_sum, 'r*', 'MarkerSize', 15, 'LineWidth', 2);
    
    title('Grid Search: Somma dei Fit vs Frequenze di Taglio');
    xlabel('f_{min} (Hz)');
    ylabel('f_{max} (Hz)');
    zlabel('Sum of Fits (%)');
    colorbar;
    view(-45, 35);
    grid on;

    % ============================================================
    % RIPRISTINO WARNING
    % ============================================================
    warning(orig_warning_state);

end