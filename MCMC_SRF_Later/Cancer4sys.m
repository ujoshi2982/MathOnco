%============================== Cancer4sys.m ==============================%
% Core 7D ODE model: dydt = f(t, y, theta)
%
% Estimated via MCMC:
%   theta = [alpha_A; alpha_B; alpha_C; alpha_D; K1; K2; K3; beta_3]
%
% Fixed (literature; edit values below):
%   Kf1, Kr1, Kf2, Kr2
%   gamma_A, gamma_AB, gamma_B, gamma_BC, gamma_C, gamma_D, gamma_I
%=========================================================================%

function dydt = Cancer4sys(~, y, theta)

    theta = theta(:);

    % ---------------- Fixed parameters (LITERATURE) ----------------
    Kf1    = 4.8;
    Kr1    = 1.8;
    Kf2    = 4.8;
    Kr2    = 1.8;

    gamma_A  = 2.1e-3;
    gamma_AB = 2.1e-3;
    gamma_B  = 2.1e-3;
    gamma_BC = 1.1e-3;
    gamma_C  = 1.1e-3;
    gamma_D  = 1.0e-3;
    gamma_I  = 1.0e-3;

    % ---------------- Estimated parameters (MCMC) -------------------
    alpha_A = theta(1);
    alpha_B = theta(2);
    alpha_C = theta(3);
    alpha_D = theta(4);
    K1      = theta(5);
    K2      = theta(6);
    K3      = theta(7);
    beta_3  = theta(8);

    % ---------------- States ----------------
    A  = y(1);
    AB = y(2);
    B  = y(3);
    BC = y(4);
    C  = y(5);
    D  = y(6);
    I  = y(7);

    % ---------------- ODEs (UNCHANGED) ----------------
    dA  = alpha_A - Kf1*A*B + (Kr1*(AB*D))/(K1 + AB) - gamma_A*A;
    dAB = Kf1*A*B - (Kr1*(AB*D))/(K1 + AB) - gamma_AB*AB;

    dB  = alpha_B - Kf1*A*B + (Kr1*(AB*D))/(K1 + AB) ...
                    - Kf2*B*C + (Kr2*(BC*AB))/(K2 + BC) - gamma_B*B;

    dBC = Kf2*B*C - (Kr2*(BC*AB))/(K2 + BC) - gamma_BC*BC;

    dC  = alpha_C - Kf2*B*C + (Kr2*(BC*AB))/(K2 + BC) - gamma_C*C;
    dD  = alpha_D - gamma_D*D;

    dI  = beta_3*(BC/(K3 + BC)) - gamma_I*I;

    dydt = [dA; dAB; dB; dBC; dC; dD; dI];

end
