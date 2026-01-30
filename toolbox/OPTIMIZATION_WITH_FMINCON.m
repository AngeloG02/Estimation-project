%% SOLUZIONE CON FMINCON - Output Error Identification
% Usa fmincon di MATLAB (robusto e testato) invece del custom optimizer



fprintf('\n=== OUTPUT ERROR IDENTIFICATION WITH FMINCON ===\n\n');

% Supponi che hai già:
% - G: function handle del modello: G(theta, freq_vect) -> [K x 2]
% - G_m: FRF misurate [K x 2] complex
% - freq_vect: vettore di frequenze [K x 1]
% - theta_0: guess iniziale [6 x 1]

%% STEP 1: Definisci la costfunction
fprintf('Setting up optimization...\n');

% Wrapper per fmincon
objective = @(theta) costFunc_outputerror(theta, G_m, freq_vect, G);

% Opzioni
options = optimoptions('fmincon', ...
    'Display', 'iter-detailed', ...          % Verbose output
    'Algorithm', 'sqp', ...                   % Sequential Quadratic Programming (robusto)
    'TolFun', 1e-8, ...                       % Tolleranza sulla costfunction
    'TolX', 1e-8, ...                         % Tolleranza sui parametri
    'MaxIterations', 200, ...
    'MaxFunctionEvaluations', 3000, ...
    'SpecifyObjectiveGradient', false, ...    % MATLAB calcola il gradiente numericamente
    'FiniteDifferenceType', 'central');       % Differenza centrale (più accurato)

%% STEP 2: Esegui l'ottimizzazione
fprintf('\nStarting optimization...\n\n');

[theta_hat, J_hat, exitflag, output, ~, grad_final, hessian_final] = fmincon(objective, theta_0, ...
    [], [], [], [], [], [], [], options);

%% STEP 3: Risultati
fprintf('\n\n========================================\n');
fprintf('OPTIMIZATION RESULTS\n');
fprintf('========================================\n\n');

fprintf('Initial cost: J(theta_0) = %.6e\n', costFunc_outputerror(theta_0, G_m, freq_vect, G));
fprintf('Final cost:   J(theta*) = %.6e\n', J_hat);
fprintf('Improvement: %.2fx\n\n', costFunc_outputerror(theta_0, G_m, freq_vect, G) / J_hat);

fprintf('Exit flag: %d\n', exitflag);
fprintf('  (1 = converged to local minimum)\n');
fprintf('  (2 = change in x too small)\n');
fprintf('  (other = check MATLAB help)\n\n');

fprintf('Iterations: %d\n', output.iterations);
fprintf('Function evaluations: %d\n', output.funcCount);
fprintf('Algorithm: %s\n\n', output.algorithm);

fprintf('Final parameters:\n');
fprintf('  theta = [');
fprintf(' %.6e', theta_hat);
fprintf(' ]\n\n');

fprintf('Gradient norm at solution: %.6e\n', norm(grad_final));
fprintf('(Should be < 1e-6 for good convergence)\n\n');

%% STEP 4: Valuta il modello ai parametri finali
fprintf('========================================\n');
fprintf('MODEL VALIDATION\n');
fprintf('========================================\n\n');

G_model = G(theta_hat, freq_vect);
H_q_model = G_model(:, 1);
H_ax_model = G_model(:, 2);

H_q_meas = G_m(:, 1);
H_ax_meas = G_m(:, 2);

err_q = H_q_meas - H_q_model;
err_ax = H_ax_meas - H_ax_model;

fprintf('H_q error:\n');
fprintf('  Absolute RMS: %.6e\n', sqrt(mean(abs(err_q).^2)));
fprintf('  Relative RMS: %.2f%%\n', 100*sqrt(mean((abs(err_q)./abs(H_q_meas)).^2)));

fprintf('\nH_ax error:\n');
fprintf('  Absolute RMS: %.6e\n', sqrt(mean(abs(err_ax).^2)));
fprintf('  Relative RMS: %.2f%%\n', 100*sqrt(mean((abs(err_ax)./abs(H_ax_meas)).^2)));

%% STEP 5: Plot di validazione
figure('Position', [100 100 1200 600]);

subplot(1, 2, 1);
semilogy(freq_vect, abs(H_q_meas), 'b.-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(freq_vect, abs(H_q_model), 'r--', 'LineWidth', 2);
grid on;
xlabel('Frequency [Hz]');
ylabel('Magnitude');
title('H_q: Pitch Rate');
legend('Measured', 'Model (fitted)', 'Location', 'best');

subplot(1, 2, 2);
semilogy(freq_vect, abs(H_ax_meas), 'b.-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(freq_vect, abs(H_ax_model), 'r--', 'LineWidth', 2);
grid on;
xlabel('Frequency [Hz]');
ylabel('Magnitude');
title('H_ax: Longitudinal Accel');
legend('Measured', 'Model (fitted)', 'Location', 'best');

sgtitle(sprintf('Output Error Identification Result: J = %.3e', J_hat));

print(gcf, 'fmincon_result.png', '-dpng', '-r150');
fprintf('\n[Plot saved: fmincon_result.png]\n');

%% HELPER FUNCTION: Cost Function
function J = costFunc_outputerror(theta, G_m, freq_vect, G)
    % Output Error cost function (normalized)
    
    K = length(freq_vect);
    
    % Evaluate model
    G_model = G(theta, freq_vect);  % [K x 2]
    
    % Normalize each output separately
    G_m_norms = max(abs(G_m), [], 1);
    G_m_norms(G_m_norms < 1e-10) = 1.0;
    
    G_m_norm = G_m ./ G_m_norms;
    G_model_norm = G_model ./ G_m_norms;
    
    % Stack complex FRF into real/imag
    Y_m = stackFRFData(G_m_norm);
    Y = stackFRFData(G_model_norm);
    
    % Error vector
    e = Y_m - Y;  % [4 x K]
    
    % Simple cost (L2 norm, unweighted for now)
    % If needed, add weighting matrix later
    J = 0.5 * sum(abs(e(:)).^2) / length(e(:));
    
    % Add small regularization to avoid NaN
    J = J + 1e-12;


    
end

%% HELPER FUNCTION: Stack FRF Data
function y = stackFRFData(g)
    % Convert complex FRF to real/imag stacked vector
    % Input: g [K x M] complex
    % Output: y [2M x K] real (each row is freq-dependent)
    
    K = size(g, 1);
    M = size(g, 2);
    
    if M == 1
        y = zeros(2, K);
        y(1, :) = real(g);
        y(2, :) = imag(g);
    elseif M == 2
        y = zeros(4, K);
        y(1, :) = real(g(:, 1));
        y(2, :) = imag(g(:, 1));
        y(3, :) = real(g(:, 2));
        y(4, :) = imag(g(:, 2));
    else
        error('g must be Kx1 or Kx2');
    end
end

fprintf('\n=== OPTIMIZATION COMPLETE ===\n\n');
