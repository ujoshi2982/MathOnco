%============================== Cancer4fun.m ==============================%
% ODE wrapper: integrates the 7-state system for a given theta over "time"
%=========================================================================%

function y = Cancer4fun(time, theta)

    % ---- Initial conditions (fixed) ----
    y0 = [0.070366; 0; 3.40861; 0; 0.121311; 0.113287; 0.001316];

    % ---- Enforce non-negativity on all 7 states ----
    options = odeset('NonNegative', 1:7);

    % ---- Solve stiff ODE system ----
    [~, y] = ode15s(@(t, y) Cancer4sys(t, y, theta), time, y0, options);

end

