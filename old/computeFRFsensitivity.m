% function dydtheta = computeFRFsensitivity( freq_vect, G, theta_0 )
% 
% %% COMPUTEFRFSENSITIVITY Compute sensitivity of the model FRF with respect to parameters.
% % 
% % input
% %   omega_vec   frequency grid
% %   G           model
% %   theta_0     parameter value
% % output
% %   dydtheta    sensitivity of FRF (real and imag parts) wrt variations of
% %               theta
% 
% % Author: 
% %   Simone Panza <simone.panza@polimi.it>
% % 
% % changelog
% %   v1
% %       first version
% 
% % perturb theta_0 by a small fraction of it 
% PERT = 0.1;     % perturb percentage
% 
% % K = length(omega_vec);
% nth = length(theta_0);
% 
% % evaluate FRF at theta_0
% G_model = G(theta_0,freq_vect);
% 
% Y0 = stackFRFData(G_model);
% 
% dydtheta = struct([]);
% 
% for ii = 1:nth
% 
%     % perturb theta
%     dtheta = PERT*abs(theta_0(ii));
%     pert_th = zeros(size(theta_0));
%     pert_th(ii) = dtheta;
%     theta = theta_0 + pert_th;   % perturbed parameter
% 
%     % evaluate FRF at theta
%     G_model = G(theta,freq_vect);
%     Y = stackFRFData(G_model);
% 
%     % compute sensitivity
%     dydtheta(ii).dy = (Y - Y0)/dtheta;
% 
% end
% 
% end

function dydtheta = computeFRFsensitivity(freq_vect, G, theta_0)
    PERT = 1e-5;  % piccolo per FD accurato
    nth = length(theta_0);
    G_model = G(theta_0, freq_vect);
    Y0 = stackFRFData(G_model);
    dydtheta = struct([]);
    
    for ii = 1:nth
        dtheta = PERT * max(abs(theta_0(ii)), 1e-3);  % avoid 0 perturbation
        
        % Central difference (più accurato)
        theta_plus = theta_0;
        theta_plus(ii) = theta_0(ii) + dtheta;
        G_plus = G(theta_plus, freq_vect);
        Y_plus = stackFRFData(G_plus);
        
        theta_minus = theta_0;
        theta_minus(ii) = theta_0(ii) - dtheta;
        G_minus = G(theta_minus, freq_vect);
        Y_minus = stackFRFData(G_minus);
        
        dydtheta(ii).dy = (Y_plus - Y_minus) / (2*dtheta);
    end
end
