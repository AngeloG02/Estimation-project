function result = optimize_fband(u, q, ax, Ts, gamma_sq_q, gamma_sq_ax, f_coh, init_sys, opt)
    % ============================================================
    % DISABILITAZIONE WARNING
    % ============================================================
    orig_warning_state = warning('query', 'all');
    warning('off', 'all');
    
    % ============================================================
    % Parametri di configurazione
    % ============================================================
    coh_thr = 0.5;
    step_Hz = 0.2;
    min_points = 8;     
    search_range = 3.0; 
    
    % ============================================================
    % 1) Banda iniziale
    % ============================================================
    valid = (gamma_sq_q > coh_thr) & (gamma_sq_ax > coh_thr);
    if ~any(valid)
        warning(orig_warning_state); 
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
    % 2) Pre-calcolo FFT su tutti i segnali 
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
    % 3) Costruzione e PRE-FILTRAGGIO della griglia (OTTIMIZZAZIONE)
    % ============================================================
    f_min_vec = max(0.2, f_min_0 - search_range) : step_Hz : f_min_0;
    f_max_vec = f_max_0 : step_Hz : min(f_nyq, f_max_0 + search_range);
    [F_MIN, F_MAX] = meshgrid(f_min_vec, f_max_vec);
    f_min_raw = F_MIN(:);
    f_max_raw = F_MAX(:);
    num_raw = length(f_min_raw);
    
    % Maschera per le combinazioni valide (f_min < f_max)
    valid_mask = (f_min_raw < f_max_raw);
    
    % Controllo dei punti minimi (fuori dal parfor)
    for k = 1:num_raw
        if valid_mask(k)
            if sum((freq_Hz >= f_min_raw(k)) & (freq_Hz <= f_max_raw(k))) < min_points
                valid_mask(k) = false;
            end
        end
    end
    
    % Estrazione delle sole combinazioni strettamente valide
    valid_indices = find(valid_mask);
    f_min_par = f_min_raw(valid_indices);
    f_max_par = f_max_raw(valid_indices);
    num_combinations = length(f_min_par);
    
    % Array pre-allocati per il ciclo parfor
    score_par  = -Inf(num_combinations, 1);
    fit_q_par  = -Inf(num_combinations, 1);
    fit_ax_par = -Inf(num_combinations, 1);
    
    fprintf('Inizio calcolo a griglia: %d combinazioni saltate, %d in elaborazione parfor...\n', ...
        num_raw - num_combinations, num_combinations);
    
    % ============================================================
    % SETUP: DataQueue per mostrare l'avanzamento
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
    % 4) Ciclo Parallelo (Grid Search sulle sole iterazioni utili)
    % ============================================================
    opt_parfor = opt;
    opt_parfor.Display = 'off';
    
    parfor i = 1:num_combinations
        warning('off', 'all'); 
        
        f_min_curr = f_min_par(i);
        f_max_curr = f_max_par(i);
        
        % Estrazione dei dati in banda
        idx = (freq_Hz >= f_min_curr) & (freq_Hz <= f_max_curr);
        U_f = U_pos(idx);
        Y_f = Y_pos(idx, :);
        W_f = W_pos(idx);
        
        % SNELLIMENTO: Creazione iddata essenziale senza metadati (più veloce)
        data = iddata(Y_f, U_f, 0, 'Frequency', W_f);
        
        try
            sys_est = greyest(data, init_sys, opt_parfor);
            
            fit_q  = sys_est.Report.Fit.FitPercent(1);
            fit_ax = sys_est.Report.Fit.FitPercent(2);
            
            % 1. Calcolo media coerenze
            idx_coh = (f_coh >= f_min_curr) & (f_coh <= f_max_curr);
            if any(idx_coh)
                mean_coh_q  = mean(gamma_sq_q(idx_coh));
                mean_coh_ax = mean(gamma_sq_ax(idx_coh));
            else
                mean_coh_q = eps; 
                mean_coh_ax = eps;
            end
            
            % 2. Media pesata dei fitting
            sum_coh = mean_coh_q + mean_coh_ax;
            if sum_coh > 0
                weighted_fit = (fit_q * mean_coh_q + fit_ax * mean_coh_ax) / sum_coh;
            else
                weighted_fit = eps;
            end
            
            % 3. Standard Deviation Relativa %
            theta_est = getpvec(sys_est, 'free');
            cov_mat   = getcov(sys_est);
            
            if isempty(cov_mat) || any(isnan(cov_mat(:)))
                mean_rel_std = Inf; 
            else
                std_theta = sqrt(diag(cov_mat));
                theta_safe = theta_est;
                theta_safe(theta_safe == 0) = eps;
                relative_std = 100 * std_theta ./ abs(theta_safe);
                mean_rel_std = mean(relative_std);
            end
            
            % 4. Assegnazione Score
            mean_rel_std = max(mean_rel_std, eps); 
            if ~isinf(mean_rel_std) && ~isnan(mean_rel_std)
                score_par(i) = weighted_fit / mean_rel_std;
            end
            
            fit_q_par(i)  = fit_q;
            fit_ax_par(i) = fit_ax;
            
        catch
            % Se fallisce, lascia i valori a -Inf
        end
        
        send(dq, i); 
    end
    
    % ============================================================
    % Ricostruzione della griglia originale per poter cercare e plottare
    % ============================================================
    score_flat = -Inf(num_raw, 1);
    score_flat(valid_indices) = score_par;
    
    fit_q_flat = -Inf(num_raw, 1);
    fit_q_flat(valid_indices) = fit_q_par;
    
    fit_ax_flat = -Inf(num_raw, 1);
    fit_ax_flat(valid_indices) = fit_ax_par;

    % ============================================================
    % 5) Individuazione del massimo globale
    % ============================================================
    [best_score, best_idx] = max(score_flat);
    if best_score == -Inf
        warning(orig_warning_state); 
        error('Nessuna combinazione ha prodotto una stima valida.');
    end
    best_f_min = f_min_raw(best_idx);
    best_f_max = f_max_raw(best_idx);
    
    % % ============================================================
    % % 6) Calcolo del modello finale (QUI reinseriamo i metadati)
    % % ============================================================
    % idx_best = (freq_Hz >= best_f_min) & (freq_Hz <= best_f_max);
    % data_best = iddata(Y_pos(idx_best, :), U_pos(idx_best), 0, 'Frequency', W_pos(idx_best));
    % data_best.InputName = 'Pitch Moment';
    % data_best.OutputName = {'Pitch rate', 'Longitudinal acceleration'};
    % data_best.OutputUnit = {'rad/s', 'm/s^2'};
    % data_best.TimeUnit = 's';
    % 
    % sys_best = greyest(data_best, init_sys, opt);
    % 
    result.f_min = best_f_min;
    result.f_max = best_f_max;
    % result.sys = sys_best;
    % result.fit_q = fit_q_flat(best_idx);
    % result.fit_ax = fit_ax_flat(best_idx);
    % result.score = best_score;
    
    % fprintf('\n========= BANDA OTTIMA TROVATA =========\n');
    % fprintf('f_min   = %.2f Hz\n', result.f_min);
    % fprintf('f_max   = %.2f Hz\n', result.f_max);
    % fprintf('fit_q   = %.2f %%\n', result.fit_q);
    % fprintf('fit_ax  = %.2f %%\n', result.fit_ax);
    % fprintf('Score   = %.4f \n', result.score);
    
    % ============================================================
    % 7) Plot 3D dei risultati
    % ============================================================
    SCORE_MATRIX = reshape(score_flat, size(F_MIN));
    SCORE_MATRIX(SCORE_MATRIX == -Inf) = NaN; 
    
    figure('Name', 'Ottimizzazione f_min e f_max', 'NumberTitle', 'off');
    surf(F_MIN, F_MAX, SCORE_MATRIX, 'EdgeColor', 'interp', 'FaceColor', 'interp');
    hold on;
    plot3(best_f_min, best_f_max, best_score, 'r*', 'MarkerSize', 15, 'LineWidth', 2);
    
    title('Grid Search: Score (Weighted Fit / % STD) vs Frequenze');
    xlabel('f_{min} (Hz)');
    ylabel('f_{max} (Hz)');
    zlabel('Score Ottimizzazione');
    colorbar;
    view(-45, 35);
    grid on;
    
    warning(orig_warning_state);
end