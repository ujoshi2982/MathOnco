%============================== Cancer4ss.m ===============================%
% Sum-of-squares objective for MCMC: compares model I(t) (state 7) vs data
%=========================================================================%

function ss = Cancer4ss(theta, data)

    theta = theta(:);
    t     = data.xdata;
    yobs  = data.ydata;

    ysim  = Cancer4fun(t, theta);
    I_sim = ysim(:, 7);

    ss = sum((I_sim - yobs).^2);

end
