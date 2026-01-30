function y = stackFRFData(g)
%% STACKFRFDATA Stack complex FRF points into real and imaginary parts.
%
% input
%   g: [Kx1] or [Kx2] complex vector/matrix of FRF data
% output
%   y: 
%       - if g is Kx1: [2xK] matrix [real(g); imag(g)]
%       - if g is Kx2: [4xK] matrix [Re(col1); Im(col1); Re(col2); Im(col2)]

K = size(g, 1);
M = size(g, 2);

if M == 1
    y = zeros(2, K);
    y(1,:) = real(g);
    y(2,:) = imag(g);
elseif M == 2
    y = zeros(4, K);
    y(1,:) = real(g(:,1));
    y(2,:) = imag(g(:,1));
    y(3,:) = real(g(:,2));
    y(4,:) = imag(g(:,2));
else
    error('g must be Kx1 or Kx2 complex matrix.');
end

end

