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
    step_Hz = 0.1;
    min_points = 8;     
    search_range = 0.7; 
    
    % ============================================================
    % 1) Banda iniziale: frequenze dove entrambe le coherence > 0.5
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
    % 2) Pre-calcolo FFT e PESI DI COERENZA SULLA GRIGLIA FFT
    % ============================================================
    N = length(u);
    fs = 1/Ts;
    f_nyq = fs/2;
    U = fft(u);
    Y = fft([q, ax]);
    Npos = floor(N/2) + 1;
    U_pos = U(1:Npos);
    Y_pos = Y(1:Npos, :);
    
    % Vettore frequenze della FFT
    freq_Hz = (0:Npos-1)' * fs / N;
    W_pos = 2*pi*freq_Hz;
    
    % -- Interpolazione della coerenza sulle frequenze della FFT --
    gamma_q_interp  = interp1(f_coh, gamma_sq_q, freq_Hz, 'linear', 'extrap');
    gamma_ax_interp = interp1(f_coh, gamma_sq_ax, freq_Hz, 'linear', 'extrap');
    
    % Saturazione fisica [0, 1] per sicurezza
    gamma_q_interp(gamma_q_interp < 0) = 0; gamma_q_interp(gamma_q_interp > 1) = 1;
    gamma_ax_interp(gamma_ax_interp < 0) = 0; gamma_ax_interp(gamma_ax_interp > 1) = 1;
    
    % Vettore dei pesi globale per greyest: W(f) = gamma^2 / (1 - gamma^2)
    gamma_mean_global = (gamma_q_interp + gamma_ax_interp) / 2;
    W_opt_global = gamma_mean_global ./ max(1 - gamma_mean_global, 0.05);
    W_opt_global = W_opt_global(:); % Assicura vettore colonna
    
    % ============================================================
    % 3) Costruzione della griglia di ricerca
    % ============================================================
    f_min_vec = max(0.01, f_min_0 - search_range) : step_Hz : f_min_0;
    f_max_vec = f_max_0 : step_Hz : min(f_nyq, f_max_0 + search_range);
    [F_MIN, F_MAX] = meshgrid(f_min_vec, f_max_vec);
    f_min_flat = F_MIN(:);
    f_max_flat = F_MAX(:);
    
    num_combinations = length(f_min_flat);
    weighted_fit_flat = -Inf(num_combinations, 1); % NUOVA METRICA
    fit_q_flat        = -Inf(num_combinations, 1);
    fit_ax_flat       = -Inf(num_combinations, 1);
    sum_var_flat     = Inf(num_combinations, 1);  
    score_flat        = -Inf(num_combinations, 1); 
    
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
    % Forziamo le opzioni di greyest identiche a quelle della stima finale
    opt_parfor = opt;
    opt_parfor.InitialState = 'backcast';
    opt_parfor.Focus = 'simulation';
    opt_parfor.EnforceStability = true;
    opt_parfor.SearchMethod = 'auto';
    opt_parfor.SearchOptions.MaxIterations = 100;
    opt_parfor.Display = 'off'; 
    
    parfor i = 1:num_combinations
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
            
            % Estraiamo i pesi esatti per la stima
            W_opt_slice = W_opt_global(idx);
            opt_temp = opt_parfor;
            opt_temp.WeightingFilter = W_opt_slice;
            
            % Estraiamo le coerenze medie del canale per questa specifica banda
            % Questo ci serve per calcolare il nuovo Score richiesto
            mean_coh_q  = mean(gamma_q_interp(idx));
            mean_coh_ax = mean(gamma_ax_interp(idx));
            
            try
                % Stima del modello CON PESATURA
                sys_est = greyest(data, init_sys, opt_temp);
                
                fit_q = sys_est.Report.Fit.FitPercent(1);
                fit_ax = sys_est.Report.Fit.FitPercent(2);
                
                % MEDIA DEI FITTING PESATA PER LA COERENZA
                % Formula della media pesata: (Fit_q * W_q + Fit_ax * W_ax) / (W_q + W_ax)
                % L'aggiunta di eps evita divisioni per zero se le coerenze dovessero essere nulle
                weighted_fit = (fit_q * mean_coh_q + fit_ax * mean_coh_ax) / (mean_coh_q + mean_coh_ax + eps);
                
                cov_mat = getcov(sys_est);
                
                if isempty(cov_mat) || any(isnan(cov_mat(:)))
                    sum_var = Inf;
                    score = -Inf;
                else
                    sum_var =sum(diag(cov_mat));
                    % NUOVO SCORE: Media pesata dei fit diviso la varianza
                    score = weighted_fit / max(sum_var, eps); 
                end
                
                weighted_fit_flat(i) = weighted_fit;
                fit_q_flat(i)        = fit_q;
                fit_ax_flat(i)       = fit_ax;
                sum_var_flat(i)     = sum_var;
                score_flat(i)        = score;
                
            catch
                % Se greyest fallisce, lasciamo i valori a -Inf
            end
        end
        send(dq, i); 
    end
    
    % ============================================================
    % 5) Individuazione del massimo globale SULLO SCORE
    % ============================================================
    [best_score, best_idx] = max(score_flat);
    if best_score == -Inf
        warning(orig_warning_state); 
        error('Nessuna combinazione ha prodotto una stima valida con matrice di covarianza calcolabile.');
    end
    
    best_f_min   = f_min_flat(best_idx);
    best_f_max   = f_max_flat(best_idx);
    best_weighted_fit = weighted_fit_flat(best_idx);
    best_var     = sum_var_flat(best_idx);
    
    % % ============================================================
    % % 6) Calcolo del modello finale sulla combinazione ottima
    % % ============================================================
    % idx_best = (freq_Hz >= best_f_min) & (freq_Hz <= best_f_max);
    % data_best = iddata(Y_pos(idx_best, :), U_pos(idx_best), 0, 'Frequency', W_pos(idx_best));
    % data_best.InputName = 'Pitch Moment';
    % data_best.OutputName = {'Pitch rate', 'Longitudinal acceleration'};
    % data_best.OutputUnit = {'rad/s', 'm/s^2'};
    % data_best.TimeUnit = 's';
    % 
    % opt_final = opt;
    % opt_final.InitialState = 'backcast';
    % opt_final.Focus = 'simulation';
    % opt_final.EnforceStability = true;
    % opt_final.SearchMethod = 'auto';
    % opt_final.Display = 'on'; 
    % opt_final.WeightingFilter = W_opt_global(idx_best);
    % 
    % fprintf('\nCalcolo del modello finale sulla banda ottima...\n');
    % sys_best = greyest(data_best, init_sys, opt_final);
    % 
    result.f_min = best_f_min;
    result.f_max = best_f_max;
    % result.sys = sys_best;
    % result.fit_q = fit_q_flat(best_idx);
    % result.fit_ax = fit_ax_flat(best_idx);
    % result.weighted_fit = best_weighted_fit;
    % result.mean_var = best_var;
    % result.score = best_score;
    % 
    % fprintf('\n========= BANDA OTTIMA TROVATA =========\n');
    % fprintf('Criterio: max( MediaFitPesata / MediaVarianze )\n');
    % fprintf('f_min             = %.2f Hz\n', result.f_min);
    % fprintf('f_max             = %.2f Hz\n', result.f_max);
    % fprintf('fit_q             = %.2f %%\n', result.fit_q);
    % fprintf('fit_ax            = %.2f %%\n', result.fit_ax);
    % fprintf('Media Fit Pesata  = %.2f %%\n', result.weighted_fit);
    % fprintf('Media varianze    = %.4e\n', result.mean_var);
    % fprintf('Score (Max)       = %.2f\n', result.score);
    % 
    % ============================================================
    % 7) Plot 3D dei risultati (Visualizza lo SCORE)
    % ============================================================
    SCORE_MATRIX = reshape(score_flat, size(F_MIN));
    SCORE_MATRIX(SCORE_MATRIX == -Inf) = NaN; 
    
    figure('Name', 'Ottimizzazione f_min e f_max', 'NumberTitle', 'off');
    surf(F_MIN, F_MAX, SCORE_MATRIX, 'EdgeColor', 'interp', 'FaceColor', 'interp');
    hold on;
    plot3(best_f_min, best_f_max, best_score, 'r*', 'MarkerSize', 15, 'LineWidth', 2);
    
    title('Grid Search: Score (MediaFitPesata / Varianza)');
    xlabel('f_{min} (Hz)');
    ylabel('f_{max} (Hz)');
    zlabel('Score Objective');
    colorbar;
    view(-45, 35);
    grid on;
    
    warning(orig_warning_state);
end