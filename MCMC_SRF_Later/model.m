%================================ model.m =================================%
% Optional: same RHS as Cancer4sys, but with explicit parameter inputs.
% (Kept as-is; useful if you want to bypass theta-packing later.)
%===========================================================================

function dydt = model(~, y, alpha_A, alpha_B, alpha_C, alpha_D, ...
                      gamma_A, gamma_B, gamma_C, gamma_D, gamma_I, gamma_AB, gamma_BC, ...
                      Kf1, Kr1, Kf2, Kr2, K1, K2, K3, beta_3)

    % [1] States
    A  = y(1);  AB = y(2);  B  = y(3);  BC = y(4);
    C  = y(5);  D  = y(6);  I  = y(7);

    % [2] ODEs (UNCHANGED)
    dA  = alpha_A - Kf1*A*B + (Kr1*(AB*D))/(K1 + AB) - gamma_A*A;
    dAB = Kf1*A*B - (Kr1*(AB*D))/(K1 + AB) - gamma_AB*AB;

    dB  = alpha_B - Kf1*A*B + (Kr1*(AB*D))/(K1 + AB) ...
                    - Kf2*B*C + (Kr2*(BC*AB))/(K2 + BC) - gamma_B*B;

    dBC = Kf2*B*C - (Kr2*(BC*AB))/(K2 + BC) - gamma_BC*BC;

    dC  = alpha_C - Kf2*B*C + (Kr2*(BC*AB))/(K2 + BC) - gamma_C*C;
    dD  = alpha_D - gamma_D*D;

    dI  = beta_3*(BC/(K3 + BC)) - gamma_I*I;

    % [3] Output
    dydt = [dA; dAB; dB; dBC; dC; dD; dI];

end


