%=========================== main_MCMC_run.m ==============================%
% Full MCMC workflow (kept identical): directory -> diary -> data -> MCMC ->
% plots -> save figures -> save workspace
%=========================================================================%

close all; clear; clc;

%% ---------------------------- Run folder/log ----------------------------
init_time = string(datetime('now','Format','ddMMyyyy_HHmmss'));
mcmc_dir  = fullfile('/mnt/data/MathOnco/MCMC_data', strcat('Model1_', init_time));
mkdir(mcmc_dir);

diary(fullfile(mcmc_dir, 'CommandWindow_Output.txt'));
disp(['Run started at ', char(init_time)]);

%% ------------------------------- Data -----------------------------------
data.xdata = [0.0000; 0.2478; 0.4625; 0.6112; 0.7433; 0.8755; 0.9746; 1.1728; 1.3214; 1.5031;
              1.6683; 1.7509; 1.8665; 2.0152; 2.1804; 2.3290; 2.5438; 2.7750; 2.9898; 3.1549;
              3.3036; 3.4523; 3.6505; 3.8322; 4.0304; 4.2451; 4.4929; 4.7737; 5.0545; 5.3518;
              5.7152; 5.9960; 6.3264; 6.5907; 6.9045; 7.2679; 7.7304; 8.0938; 8.5398; 9.0849;
              9.6630; 10.1255; 10.6376; 11.2322; 11.7939; 12.4711; 12.9666; 13.4952; 14.2055; 15.1635;
              16.0059; 16.4519; 16.6667];

data.ydata = [0.0000; 0.0153; 0.1070; 0.2190; 0.3616; 0.5348; 0.6570; 0.9728; 1.2275; 1.5535;
              1.8489; 2.0221; 2.2054; 2.4703; 2.7250; 2.9491; 3.2649; 3.5195; 3.7538; 3.9066;
              4.0390; 4.1613; 4.2835; 4.3854; 4.4771; 4.5688; 4.6299; 4.7114; 4.7725; 4.8234;
              4.8744; 4.8947; 4.9253; 4.9355; 4.9559; 4.9559; 4.9660; 4.9864; 4.9864; 4.9966;
              5.0068; 5.0068; 5.0068; 5.0068; 5.0068; 5.0068; 5.0170; 5.0272; 5.0272; 5.0170;
              5.0272; 5.0272; 5.0170];

%% ----------------------------- MCMC Setup -------------------------------
model.ssfun = @Cancer4ss;
model.N     = length(data.ydata);

Kini = [0.01; 1.12; 0.02; 0.02; 1.5; 1.9; 1.5; 1.2; 0.01; 0.01];
K0   = Kini;

params = {
    {'alpha_A',  K0(1),  1e-6, 1e1, K0(1),  1}
    {'alpha_B',  K0(2),  1e-6, 1e1, K0(2),  1}
    {'alpha_C',  K0(3),  1e-6, 1e1, K0(3),  1}
    {'alpha_D',  K0(4),  1e-6, 1e1, K0(4),  1}
    {'K1',       K0(5),  1e-6, 1e2, K0(5),  1}
    {'K2',       K0(6),  1e-6, 1e2, K0(6),  1}
    {'K3',       K0(7),  1e-6, 1e2, K0(7),  1}
    {'beta_3',   K0(8),  1e-6, 1e2, K0(8),  1}
%     {'Kr1',      K0(9),  1e-6, 1e2, K0(9),  1}
%     {'Kr2',      K0(10), 1e-6, 1e2, K0(10), 1}
    {'gamma_AB', K0(9), 1e-6, 1e2, K0(9), 1}
    {'gamma_BC', K0(10), 1e-6, 1e2, K0(10), 1}
};

%% ------------------------------ MCMC Run --------------------------------
options = struct('nsimu', 5e4, 'updatesigma', 4, 'verbosity', 1, 'waitbar', 1, 'burnintime', 5e3);
[results, chain, s2chain] = mcmcrun(model, data, params, options);

options.nsimu     = 5e5;
options.burnintime = 5e4;
[results, chain, s2chain] = mcmcrun(model, data, params, options, results);

%% ------------------------------ Diagnostics -----------------------------
figure(1); clf;
mcmcplot(chain, [], results, 'hist');

figure(2); clf;
mcmcplot(chain, [], results, 'denspanel', 2);

figure(3); clf;
mcmcplot(chain, [], results);
chainstats(chain, results);

%% ----------------------- Model fit using mean(theta) --------------------
y0          = [0.070366; 0; 3.40861; 0; 0.121311; 0.113287; 0.001316];
ode_options = odeset('NonNegative', 1:7);

[t, y] = ode15s(@(t, y) Cancer4sys(t, y, mean(chain)), data.xdata, y0, ode_options);

figure(4); clf;
plot(t, y(:,7), 'b-', 'LineWidth', 1.5, 'DisplayName', 'PI3K (Simulated)');
hold on;
scatter(data.xdata, data.ydata, 50, 'r', 'filled', 'DisplayName', 'Observed Data');
hold off;

legend('Location', 'Best');
xlabel('Time (Minutes)');
ylabel('Concentration (\muM)');
title('Simulated vs Observed Data for PI3K');
grid on;

%% ---------------------------- Save figures ------------------------------
figureNames = {
    'MCMC_Histograms'
    'MCMC_DensityPanel'
    'MCMC_ChainTraces'
    'Simulated_vs_Observed_PI3K'
};

figHandles = findall(0, 'Type', 'figure');

for i = 1:min(length(figHandles), length(figureNames))
    fig      = figHandles(i);
    fig_name = figureNames{i};
    saveas(fig, fullfile(mcmc_dir, strcat(fig_name, '.fig')));
    saveas(fig, fullfile(mcmc_dir, strcat(fig_name, '.svg')));
end

%% ------------------------ Save workspace and exit -----------------------
close all;
save(fullfile(mcmc_dir, sprintf('MCMC_Workspace_%s.mat', init_time)), '-v7.3');

disp('MCMC run complete.');
diary off;

