function H = Model_system_FRF(theta,freq_vect)
% restituisce le due FRF (uscita1 = q, uscita2 = a_x)

    theta = theta(:);
    X_u = theta(1);
    X_q = theta(2);
    M_u = theta(3);
    M_q = theta(4);
    X_d = theta(5);
    M_d = theta(6);

    g = 9.81;

    A = [X_u, X_q, -g;
         M_u, M_q,  0;
           0,   1,  0];

    B = [X_d;
         M_d;
           0];

    C = [0,   1, 0;      % prima uscita (q)
         X_u, X_q, 0];   % seconda uscita (a_x)

    D = [0; X_d];        % [2x1] feedthrough

    K = numel(freq_vect);
    H_q  = zeros(K,1);
    H_ax = zeros(K,1);

    I3 = eye(3);
    for k = 1:K
        s  = 1j*freq_vect(k);
        Gs = C * ((s*I3 - A) \ B) + D;   % [2x1]
        H_q(k)  = Gs(1);
        H_ax(k) = Gs(2);
    end

    H =  [H_q, H_ax];
end

