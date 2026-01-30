function f = computeFit( g_m, omega_vec, G)
%% COMPUTEFIT Compute fit function of model vs measured FRF.
% 
% input
%   g_m         measured (nonparametric) FRF data
%   omega_vec   frequency grid [rad/s]
%   G           tf model
% 
% output
%   f           fit of identified model FRF to the nonparametric FRF

% Author: 
%   Simone Panza <simone.panza@polimi.it>
% 
% changelog
%   v1
%       first version

% FRF
g = squeeze(freqresp(G, omega_vec));  

% error
e = g_m - g;

% compute FIT 
f = 1 - sum(abs(e).^2) / sum(abs(g_m).^2);


end

