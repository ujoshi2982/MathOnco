%============================== Cancer4sys.m ==============================%
% Core 7D ODE model: dydt = f(t,y,theta)
% theta = [alpha_A, alpha_B, alpha_C, alpha_D, K1, K2, K3, beta_3,
%          Kr1, Kr2, gamma_AB, gamma_BC]
%=========================================================================%

function dydt = Cancer4sys(~, y, theta)

    % ---- Fixed decay parameters ----
    gamma_A = 2.1e-3;
    gamma_B = 2.1e-3;
    gamma_C = 1.1e-3;
    gamma_D = 1.0e-3;
    gamma_I = 1.0e-3;
%     gamma_AB = 2.1e-3;
%     gamma_BC = 1.1e-3;
    

    % ---- Fixed forward rates ----
    Kf1 = 4.8;
    Kf2 = 4.8;
    Kr1 = 0.6;
    Kr2 = 0.6;

    % ---- Estimated parameters (theta) ----
    alpha_A  = theta(1);
    alpha_B  = theta(2);
    alpha_C  = theta(3);
    alpha_D  = theta(4);
    K1       = theta(5);
    K2       = theta(6);
    K3       = theta(7);
    beta_3   = theta(8);
%     Kr1      = theta(9);
%     Kr2      = theta(10);
    gamma_AB = theta(9);
    gamma_BC = theta(10);

    % ---- States ----
    A  = y(1);
    AB = y(2);
    B  = y(3);
    BC = y(4);
    C  = y(5);
    D  = y(6);
    I  = y(7);

    % ---- ODEs ----
    dA_dt  = alpha_A - Kf1*A*B + (Kr1*(AB*D))/(K1 + AB) - gamma_A*A;

    dAB_dt = Kf1*A*B - (Kr1*(AB*D))/(K1 + AB) - gamma_AB*AB;

    dB_dt  = alpha_B - Kf1*A*B + (Kr1*(AB*D))/(K1 + AB) - Kf2*B*C + (Kr2*(BC*AB))/(K2 + BC) - gamma_B*B;

    dBC_dt = Kf2*B*C - (Kr2*(BC*AB))/(K2 + BC) - gamma_BC*BC;

    dC_dt  = alpha_C - Kf2*B*C + (Kr2*(BC*AB))/(K2 + BC) - gamma_C*C;

    dD_dt  = alpha_D - gamma_D*D;

    dI_dt  = beta_3 * (BC/(K3 + BC)) - gamma_I*I; 

    dydt = [dA_dt; dAB_dt; dB_dt; dBC_dt; dC_dt; dD_dt; dI_dt];

end

