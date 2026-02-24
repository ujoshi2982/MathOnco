%============================== Cancer4fun.m ==============================%
% ODE wrapper (7-state): integrates Cancer4sys over given time for theta
%=========================================================================%

function y = Cancer4fun(time, theta)

    theta = theta(:);

    y0 = [0.07035; 0; 3.40861; 0; 0.12131; 0.11329; 0.001316];

    ode_opts = odeset('NonNegative', 1:7, 'RelTol',1e-6, 'AbsTol',1e-9);

    [~, y] = ode15s(@(t, y) Cancer4sys(t, y, theta), time, y0, ode_opts);

end
