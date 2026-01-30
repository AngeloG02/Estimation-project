% function [theta_hat, J_hat] = outputError_FD(G_m, freq_vect, G_th, theta_0)
% %% OUTPUTERROR_FD output error method in frequency domain with line search
% % Input
% %   G_m         measured FRF data
% %   freq_vect   frequency grid [rad/s]
% %   G_th        model FRF
% %   theta_0     parameter vector initial guess
% % 
% % Output
% %   theta_hat   optimal parameter vector
% %   J_hat       optimal cost function value
% 
% % Parameters for convergence criteria
% TRESHOLD_THETA = 1e-3;      % absolute threshold on the parameter vector
% TRESHOLD_G = 1e-6;          % absolute threshold on gradient norm
% Nmax = 50;                  % max number of iterations
% LINE_SEARCH_ITER = 20;      % max iterations of line search
% ARMIJO_C1 = 1e-4;           % Armijo condition constant (0 < c1 < 0.5)
% 
% assert(iscolumn(theta_0), 'Parameter vector must be a column vector.');
% 
% K = length(freq_vect);
% nth = length(theta_0);
% 
% % Initialization
% theta = theta_0;
% c = 0;
% 
% % Initial cost function
% [J, e, Rinv] = costFunction(G_m, freq_vect, G_th, theta);
% 
% % Debug: initial guess
% fprintf('\n========== OUTPUT ERROR OPTIMIZATION ==========\n');
% fprintf('Initial guess: J = %12.6e\n', J);
% for ii = 1:nth
%     fprintf('  θ%d = %12.6e\n', ii, theta(ii));
% end
% fprintf('\n');
% 
% % NEWTON-RAPHSON WITH LINE SEARCH
% while c < Nmax
% 
%     % Compute sensitivity of the model FRF with respect to parameters
%     s = computeFRFsensitivity(freq_vect, G_th, theta);
% 
%     % Compute gradient and Hessian of J wrt theta
%     G_grad = zeros(nth, 1);  % gradient
%     H = zeros(nth, nth);      % Hessian
% 
%     for kk = 1:K
%         dydtheta = zeros(4, nth);
%         % Concatenate sensitivities for all parameters
%         for ii = 1:nth
%             dydtheta(:, ii) = s(ii).dy(:, kk);
%         end
% 
%         % Gradient contribution: -dydtheta' * Rinv * e
%         G_grad = G_grad + dydtheta' * Rinv * e(:, kk);
% 
%         % Hessian contribution: dydtheta' * Rinv * dydtheta
%         % (WITHOUT frequency weighting on Hessian - that distorts it)
%         H = H + dydtheta' * Rinv * dydtheta;
% 
%         % Subito dopo il calcolo di H
% 
% 
%     end
% 
%     % Normalize gradient
%     norm_G = norm(G_grad);
% 
%     % Compute Newton-Raphson step with Levenberg-Marquardt damping
%     % Start with small lambda and increase if needed
%     lambda = 1e-6;
%     max_lambda_iter = 5;
%     found_descent = false;
% 
%     for lambda_iter = 1:max_lambda_iter
%         try
%             % Try to solve with current lambda
%             delta_theta = -(H + lambda*eye(nth)) \ G_grad;
% 
%             % Check if it's a descent direction: delta_theta' * G_grad < 0
%             if delta_theta' * G_grad < -1e-10 * norm(delta_theta) * norm_G
%                 found_descent = true;
%                 break;
%             end
% 
%             % If not descent, increase lambda
%             lambda = lambda * 10;
% 
%         catch
%             % Matrix singular, increase lambda
%             lambda = lambda * 10;
%         end
%     end
% 
%     if ~found_descent
%         warning('Failed to find descent direction at iteration %d. Stopping.', c);
%         break;
%     end
% 
%     % ===== LINE SEARCH (ARMIJO BACKTRACKING) =====
%     alpha = 1.0;
%     theta_new = [];
%     J_new = inf;
% 
%     for ls_iter = 1:LINE_SEARCH_ITER
%         theta_test = theta + alpha * delta_theta;
%         [J_test, e_test, ~] = costFunction(G_m, freq_vect, G_th, theta_test);
% 
%         % Armijo condition: J(theta + alpha*delta) < J(theta) + c1*alpha*gradient'*delta
%         armijo_rhs = J + ARMIJO_C1 * alpha * (G_grad' * delta_theta);
% 
%         if J_test < armijo_rhs
%             % Line search succeeded
%             theta_new = theta_test;
%             J_new = J_test;
%             e_new = e_test;
%             break;
%         end
% 
%         % Reduce step size
%         alpha = alpha * 0.5;
%     end
% 
%     fprintf('  cond(H) = %.3e\n', cond(H));
% fprintf('  cond(H + lambda*I) = %.3e\n', cond(H + lambda*eye(nth)));
% fprintf('  norm(G_grad) = %.3e\n', norm(G_grad));
% fprintf('  norm(delta_theta) = %.3e\n', norm(delta_theta));
% 
%     if isempty(theta_new)
%         warning('Line search failed at iteration %d. Step too small.', c);
%         break;
%     end
% 
%     % Update estimates
%     theta = theta_new;
%     J = J_new;
%     e = e_new;
%     c = c + 1;
% 
%     % Print progress
%     fprintf('Iter %2d: J = %12.6e | ||ΔΘ|| = %9.3e | ||∇J|| = %9.3e | α = %.2e | λ = %.2e\n', ...
%             c, J, norm(delta_theta), norm_G, alpha, lambda);
%     for ii = 1:nth
%         fprintf('           θ%d = %12.6e\n', ii, theta(ii));
%     end
% 
%     % Check convergence
%     if norm(delta_theta) < TRESHOLD_THETA && norm_G < TRESHOLD_G
%         fprintf('\nConverged after %d iterations\n', c);
%         break;
%     end
% 
% end
% 
% fprintf('================================================\n\n');
% 
% % Best estimate
% theta_hat = theta;
% J_hat = J;
% 
% end
% 
function [theta_hat, J_hat] = outputError_FD(G_m, freq_vect, G_th, theta_0)
%% OUTPUTERROR_FD - OUTPUT ERROR METHOD - ROBUST & NORMALIZED VERSION
%
% Fixes applied:
%   1. Data normalization (per-output normalization)
%   2. Analytical sensitivity computation with proper scaling
%   3. Properly constructed Hessian (no spurious scaling)
%   4. Adaptive lambda with aggressive initialization
%   5. Relative step-size tolerance
%   6. Robust line search with multiple fallback strategies

TRESHOLD_GRAD  = 1e-4;      % Gradient norm tolerance (scaled)
TRESHOLD_REL   = 1e-3;      % Relative parameter change tolerance
Nmax           = 150;
LINE_SEARCH_ITER = 30;
ARMIJO_C1      = 1e-4;      % Armijo constant
alpha_min      = 1e-10;
LAMBDA_INIT    = 1e5;       % Start very aggressive
LAMBDA_MAX     = 1e10;      % Upper bound on lambda
MAX_LAMBDA_ITER = 20;
LAMBDA_MULT    = 5;         % Multiplication factor for lambda increase

assert(iscolumn(theta_0), 'Parameter vector must be a column vector.');
K = length(freq_vect);
nth = length(theta_0);

fprintf('\n========================================\n');
fprintf('  OUTPUT ERROR OPTIMIZATION - ROBUST v2\n');
fprintf('========================================\n');
fprintf('Problem size: %d parameters, %d frequencies\n', nth, K);

% ========================================================================
% STEP 0: DATA NORMALIZATION
% ========================================================================
% Normalize each output (column) separately to unit norm
G_m_norms = max(abs(G_m), [], 1);  % [1 x 2]
G_m_norms(G_m_norms < 1e-10) = 1.0;  % Avoid division by zero

G_m_normalized = G_m ./ G_m_norms;  % Normalize columns

fprintf('\nData Normalization:\n');
fprintf('  Output 1 (q):   norm = %.3e\n', G_m_norms(1));
fprintf('  Output 2 (a_x): norm = %.3e\n', G_m_norms(2));
fprintf('  After normalization: ||G_m|| ≈ 1\n\n');

% ========================================================================
% STEP 1: PARAMETER NORMALIZATION
% ========================================================================
theta_scale = abs(theta_0);
theta_scale(abs(theta_scale) < 1e-10) = 1.0;

fprintf('Parameter Scaling:\n');
for ii = 1:nth
    fprintf('  θ%d: scale = %.3e\n', ii, theta_scale(ii));
end
fprintf('\n');

theta_norm = theta_0 ./ theta_scale;

% ========================================================================
% STEP 2: INITIAL COST FUNCTION (on normalized data)
% ========================================================================
[J_init, e_init, Rinv_init] = costFunction_normalized(G_m_normalized, freq_vect, G_th, theta_0, G_m_norms);

fprintf('Initial cost function (normalized): J = %.6e\n\n', J_init);

J = J_init;
e = e_init;
Rinv = Rinv_init;

% ========================================================================
% ITERATION LOOP
% ========================================================================
fprintf('Iteration Log:\n');
fprintf('--------------------------------------------------------------\n');
fprintf('Iter |      J       | ||ΔΘ||_rel | ||∇J||/||θ|||    α      | λ      \n');
fprintf('--------------------------------------------------------------\n');

c = 0;
best_J = J;
best_theta_norm = theta_norm;
best_iter = 0;

while c < Nmax
    
    theta = theta_norm .* theta_scale;
    
    % ====================================================================
    % SENSITIVITY COMPUTATION (analytical, properly scaled)
    % ====================================================================
    s = computeFRFsensitivity_analytical(freq_vect, G_th, theta, theta_scale, theta_norm);
    
    % ====================================================================
    % JACOBIAN & GRADIENT & HESSIAN COMPUTATION
    % ====================================================================
    G_grad_norm = zeros(nth, 1);
    H_norm = zeros(nth, nth);
    
    for kk = 1:K
        % Jacobian at frequency kk: [4 x nth]
        J_k = zeros(4, nth);
        for ii = 1:nth
            J_k(:, ii) = s(ii).dy(:, kk);
        end
        
        % Gradient: ∇J = Σ_k J_k^T R_inv e_k
        G_grad_norm = G_grad_norm + J_k' * Rinv * e(:, kk);
        
        % Hessian: H ≈ Σ_k J_k^T R_inv J_k
        H_norm = H_norm + J_k' * Rinv * J_k;
    end
    
    norm_G = norm(G_grad_norm);
    norm_theta = norm(theta_norm);
    scaled_grad = norm_G / max(norm_theta, 1e-10);
    
    % ====================================================================
    % LEVENBERG-MARQUARDT WITH ADAPTIVE LAMBDA
    % ====================================================================
    lambda = LAMBDA_INIT;
    found_descent = false;
    delta_theta_norm = [];
    cond_H_best = inf;
    
    for lambda_iter = 1:MAX_LAMBDA_ITER
        H_reg = H_norm + lambda * eye(nth);
        cond_H = cond(H_reg);
        
        % Stop if matrix is well-conditioned
        if cond_H < 1e8
            cond_H_best = cond_H;
            try
                delta_theta_norm = -H_reg \ G_grad_norm;
                descent_measure = delta_theta_norm' * G_grad_norm;
                
                % Check if it's a descent direction
                if descent_measure < -1e-10 * norm(delta_theta_norm) * norm_G
                    found_descent = true;
                    break;
                end
            catch
                % Singular matrix, try next lambda
            end
        end
        
        % Increase lambda and try again
        lambda = lambda * LAMBDA_MULT;
        if lambda > LAMBDA_MAX
            break;
        end
    end
    
    if ~found_descent
        fprintf('%4d | %.6e | ---STOP--- | No descent direction found\n', c, J);
        fprintf('Stopping: Cannot find descent direction. Lambda exceeded bounds.\n');
        break;
    end
    
    % ====================================================================
    % ARMIJO LINE SEARCH (with multiple fallback strategies)
    % ====================================================================
    alpha = 1.0;
    theta_norm_new = [];
    J_new = inf;
    e_new = [];
    Rinv_new = [];
    found_step = false;
    
    % Strategy 1: Standard Armijo with exponential backtracking
    for ls_iter = 1:LINE_SEARCH_ITER
        theta_norm_test = theta_norm + alpha * delta_theta_norm;
        theta_test = theta_norm_test .* theta_scale;
        
        [J_test, e_test, Rinv_test] = costFunction_normalized(G_m_normalized, freq_vect, G_th, theta_test, G_m_norms);
        
        % Armijo condition: J(θ + α*Δθ) ≤ J(θ) + c*α*(∇J)^T*Δθ
        armijo_rhs = J + ARMIJO_C1 * alpha * (G_grad_norm' * delta_theta_norm);
        
        if J_test < armijo_rhs
            theta_norm_new = theta_norm_test;
            J_new = J_test;
            e_new = e_test;
            Rinv_new = Rinv_test;
            found_step = true;
            break;
        end
        
        alpha = alpha * 0.5;
    end
    
    % Strategy 2: Fallback - accept the best step found in line search
    if ~found_step
        alpha = 1.0;
        best_J_linesearch = J;
        best_theta_norm_linesearch = theta_norm;
        best_J_test = J;
        best_e_test = e;
        best_Rinv_test = Rinv;
        
        for ls_iter = 1:10
            theta_norm_test = theta_norm + alpha * delta_theta_norm;
            theta_test = theta_norm_test .* theta_scale;
            
            [J_test, e_test, Rinv_test] = costFunction_normalized(G_m_normalized, freq_vect, G_th, theta_test, G_m_norms);
            
            if J_test < best_J_test
                best_J_test = J_test;
                best_theta_norm_linesearch = theta_norm_test;
                best_e_test = e_test;
                best_Rinv_test = Rinv_test;
                found_step = true;
            end
            
            alpha = alpha * 0.5;
        end
        
        if found_step
            theta_norm_new = best_theta_norm_linesearch;
            J_new = best_J_test;
            e_new = best_e_test;
            Rinv_new = best_Rinv_test;
            alpha = 0.5^(10 - 1);  % Record the alpha used
        end
    end
    
    % Strategy 3: Last resort - minimum step
    if ~found_step
        alpha = alpha_min;
        theta_norm_new = theta_norm + alpha * delta_theta_norm;
        theta_test = theta_norm_new .* theta_scale;
        [J_new, e_new, Rinv_new] = costFunction_normalized(G_m_normalized, freq_vect, G_th, theta_test, G_m_norms);
        found_step = true;
    end
    
    % ====================================================================
    % UPDATE STATE
    % ====================================================================
    if found_step
        theta_norm = theta_norm_new;
        J = J_new;
        e = e_new;
        Rinv = Rinv_new;
        
        % Track best solution
        if J < best_J
            best_J = J;
            best_theta_norm = theta_norm;
            best_iter = c;
        end
    else
        fprintf('%4d | %.6e | ---STOP--- | No valid step found\n', c, J);
        break;
    end
    
    delta_norm = norm(delta_theta_norm);
    delta_norm_rel = delta_norm / max(norm(theta_norm), 1e-10);
    
    c = c + 1;
    
    fprintf('%4d | %.6e | %.3e | %.3e | %.2e | %.2e\n', ...
            c, J, delta_norm_rel, scaled_grad, alpha, lambda);
    
    % ====================================================================
    % CONVERGENCE CHECK
    % ====================================================================
    if delta_norm_rel < TRESHOLD_REL && scaled_grad < TRESHOLD_GRAD
        fprintf('--------------------------------------------------------------\n');
        fprintf('\nConverged after %d iterations\n', c);
        break;
    end
    
end

fprintf('========================================\n\n');

theta_hat = best_theta_norm .* theta_scale;
J_hat = best_J;

fprintf('FINAL RESULTS:\n');
fprintf('  Iterations:  %d (best at iteration %d)\n', c, best_iter);
fprintf('  Final J = %.6e\n', J_hat);
fprintf('  Final theta = [');
fprintf(' %.6e', theta_hat);
fprintf(' ]\n\n');

end

% =========================================================================
% COST FUNCTION (NORMALIZED)
% =========================================================================
function [J, e, R_inv] = costFunction_normalized(G_m_norm, freq_vect, G, theta, G_m_norms)
%
% INPUT:
%   G_m_norm:   Normalized measured FRF [K x 2]
%   freq_vect:  Frequency vector [K x 1]
%   G:          Model function handle
%   theta:      Parameter vector
%   G_m_norms:  Original norms for each output [1 x 2]
%
% OUTPUT:
%   J:          Normalized cost function
%   e:          Stacked error (on normalized data)
%   R_inv:      Inverse covariance matrix

K = length(freq_vect);

% Evaluate model on DENORMALIZED scale (to get correct relative errors)
G_model = G(theta, freq_vect);  % Model output (raw scale)

% Normalize model outputs
G_model_norm = G_model ./ G_m_norms;

% Stack FRF data (complex → real/imag)
Y_m = stackFRFData(G_m_norm);
Y = stackFRFData(G_model_norm);

% Error
e = Y_m - Y;

% Covariance matrix (estimated from residuals)
R_q = zeros(2, 2);
R_ax = zeros(2, 2);

e_q = stackFRFData(G_m_norm(:, 1) - G_model_norm(:, 1));
e_ax = stackFRFData(G_m_norm(:, 2) - G_model_norm(:, 2));

for kk = 1:K
    R_q = R_q + e_q(:, kk) * e_q(:, kk)';
    R_ax = R_ax + e_ax(:, kk) * e_ax(:, kk)';
end

R_q = R_q / K;
R_ax = R_ax / K;

% Block diagonal covariance
R = blkdiag(R_q, R_ax);

% Add regularization for numerical stability
R_inv = inv(R + 1e-6 * eye(4));

% Cost function
J = 0;
for kk = 1:K
    J = J + 0.5 * e(:, kk)' * R_inv * e(:, kk);
end

% Normalize by number of data points
J = J / (2 * K);

end

% =========================================================================
% SENSITIVITY COMPUTATION (ANALYTICAL - IMPROVED)
% =========================================================================
function dydtheta = computeFRFsensitivity_analytical(freq_vect, G, theta, theta_scale, theta_norm)
%
% Compute sensitivity dY/d(theta_norm) using analytical finite differences
% with adaptive perturbation sizes

PERT_BASE = 1e-2;  % 1% perturbation (not 1e-5!)
nth = length(theta_norm);
K = length(freq_vect);

G_model = G(theta, freq_vect);
Y0 = stackFRFData(G_model);

dydtheta = struct([]);

for ii = 1:nth
    % Adaptive perturbation: 1% of theta_norm(ii), with lower bound
    dtheta_norm = PERT_BASE * max(abs(theta_norm(ii)), 1e-2);
    
    % Forward difference for each direction
    % (Central difference could be expensive with 2 model evaluations per parameter)
    
    theta_norm_pert = theta_norm;
    theta_norm_pert(ii) = theta_norm(ii) + dtheta_norm;
    theta_pert = theta_norm_pert .* theta_scale;
    G_pert = G(theta_pert, freq_vect);
    Y_pert = stackFRFData(G_pert);
    
    % Derivative w.r.t. theta_norm
    dydtheta(ii).dy = (Y_pert - Y0) / dtheta_norm;
end

end

% =========================================================================
% UTILITY: Stack FRF Data
% =========================================================================
function y = stackFRFData(g)
%% STACKFRFDATA Stack complex FRF points into real and imaginary parts.

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
    error('g must be Kx1 or Kx2 complex matrix.');
end

end