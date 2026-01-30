function [ J, e,R_inv ] = costFunction( G_m, freq_vect, G, theta)
%% COST FUNCTION Evaluate cost function and error
% 
% input
%   g_m         measured FRF data
%   omega_vec   frequency grid [rad/s]
%   G           tf model parametrized in theta
%   theta       parameter vector
% 
% output
%   J           cost function value
%   e           estimation error as a function of frequency

% Author: 
%   Simone Panza <simone.panza@polimi.it>
% 
% changelog
%   v1
%       first version

% Stack complex FRF points into a vector of real and imaginary parts
Y_m = stackFRFData(G_m);

K = length(freq_vect);
% nth = length(theta);

G_model = G(theta,freq_vect); % model FRF evaluate at freq_vect

% stak the model FRF
Y = stackFRFData(G_model);

% error
e = Y_m - Y;

% cost function computation
J = 0;
% 
e_q  = stackFRFData( G_m(:,1)  - G_model(:,1));
e_ax = stackFRFData( G_m(:,2)  - G_model(:,2));
% 
% % varianze per le due uscite (stesso valore per Re/Im)
% var_q  = mean(abs(e_q ).^2) / 2;
% var_ax = mean(abs(e_ax).^2) / 2;

% R_diag = [var_q*ones(2*K,1); var_ax*ones(2*K,1)];

R_q = zeros(2,2);
R_ax = zeros(2,2);
for kk = 1:K
    R_q = R_q + e_q(:,kk)'*e_q(:, kk);
    R_ax = R_ax + e_ax(:,kk)'*e_ax(:, kk);
end
R_q = R_q/K;
R_ax = R_ax/K;


% Opzione 1: usa la matrice di covarianza completa
R = blkdiag(R_q, R_ax);  % mantieni correlazioni
R_inv = inv(R + 1e-10*eye(4));  % regularizzazione

% % Opzione 2: se vuoi matrice diagonale, fallo correttamente
% var_q = (trace(R_q))/(2*K);     % media varianze
% var_ax = (trace(R_ax))/(2*K);
% R = blkdiag(var_q*eye(2), var_ax*eye(2));
% R_inv = inv(R);

for kk = 1:K
    J = J + 1/2*e(:,kk)'*R_inv*e(:, kk);
end

end
