%============================== Cancer4ss.m ===============================%
% Sum-of-squares objective for MCMC: compares model I(t) (state 7) vs data
%=========================================================================%

function ss = Cancer4ss(theta, data)

    % ---- Experimental data ----
    time  = data.xdata;
    ydata = data.ydata;

    % ---- Model output ----
    ymodel  = Cancer4fun(time, theta);
    ymodel2 = ymodel(:, 7);      % state-7 = I

    % ---- SSE ----
    ss = sum((ymodel2 - ydata).^2);

end

