%=========================== main_MCMC_run.m ==============================%
% Full MCMC workflow:
%   [1] Create run directory + start log
%   [2] Load/define data
%   [3] Configure MCMC model + parameters
%   [4] Run MCMC (Chain 1 + Chain 2 continuation)
%   [5] Publication-style diagnostics plots (ONE figure, 2x4 panels) using mcmcplot
%   [6] Model fit using mean(theta) (post-burn)
%   [6.4] Posterior predictive band (median + 95% CI) in a separate figure
%   [7] Save figures + workspace and exit
%=========================================================================%

close all; clear; clc;

%% ---------------------------- Run folder/log ----------------------------
init_time = string(datetime('now','Format','ddMMyyyy_HHmmss'));
mcmc_dir  = fullfile('/mnt/data/MathOnco/MCMC_data', strcat('Model1_', init_time));
mkdir(mcmc_dir);

diary(fullfile(mcmc_dir, 'CommandWindow_Output.txt'));
disp(['Run started at ', char(init_time)]);

%% [2] ------------------------------- Data -------------------------------
data.xdata = [0.0000; 0.2478; 0.4625; 0.6112; 0.7433; 0.8755; 0.9746; 1.1728; 1.3214; 1.5031;
              1.6683; 1.7509; 1.8665; 2.0152; 2.1804; 2.3290; 2.5438; 2.7750; 2.9898; 3.1549;
              3.3036; 3.4523; 3.6505; 3.8322; 4.0304; 4.2451; 4.4929; 4.7737; 5.0545; 5.3518;
              5.7152; 5.9960; 6.3264; 6.5907; 6.9045; 7.2679; 7.7304; 8.0938; 8.5398; 9.0849;
              9.6630; 10.1255; 10.6376; 11.2322; 11.7939; 12.4711; 12.9666; 13.4952; 14.2055;
              15.1635; 16.0059; 16.4519; 16.6667];

data.ydata = [0.0000; 0.0153; 0.1070; 0.2190; 0.3616; 0.5348; 0.6570; 0.9728; 1.2275; 1.5535;
              1.8489; 2.0221; 2.2054; 2.4703; 2.7250; 2.9491; 3.2649; 3.5195; 3.7538; 3.9066;
              4.0390; 4.1613; 4.2835; 4.3854; 4.4771; 4.5688; 4.6299; 4.7114; 4.7725; 4.8234;
              4.8744; 4.8947; 4.9253; 4.9355; 4.9559; 4.9559; 4.9660; 4.9864; 4.9864; 4.9966;
              5.0068; 5.0068; 5.0068; 5.0068; 5.0068; 5.0068; 5.0170; 5.0272; 5.0272; 5.0170;
              5.0272; 5.0272; 5.0170];

%% [3] ---------------------------- MCMC Setup ----------------------------
model.ssfun = @Cancer4ss;
model.N     = length(data.ydata);

Kini = [0.01; 1.12; 0.02; 0.02; 1.5; 1.9; 1.5; 1.2];
K0   = Kini;

params = {
    {'alpha_A', K0(1), 1e-6, 1e1, K0(1), 1}
    {'alpha_B', K0(2), 1e-6, 1e1, K0(2), 1}
    {'alpha_C', K0(3), 1e-6, 1e1, K0(3), 1}
    {'alpha_D', K0(4), 1e-6, 1e1, K0(4), 1}
    {'K1',      K0(5), 1e-6, 1e1, K0(5), 1}
    {'K2',      K0(6), 1e-6, 1e1, K0(6), 1}
    {'K3',      K0(7), 1e-6, 1e1, K0(7), 1}
    {'beta_3',  K0(8), 1e-6, 1e1, K0(8), 1}
};

param_names = cellfun(@(c) c{1}, params, 'UniformOutput', false);

%% [4] ------------------------------ MCMC Run ----------------------------
% NOTE: warmup + continuation (not two independent chains)

nsimu1 = 5e4;  burn1 = 5e3;
nsimu2 = 5e5;  burn2 = 5e4;   % change to 5e5 / 5e6 if desired

options = struct( ...
    'nsimu',       nsimu1, ...
    'updatesigma', 4, ...
    'verbosity',   1, ...
    'waitbar',     1, ...
    'burnintime',  burn1);

fprintf('Starting Chain 1 (warmup) (%g iterations)...\n', nsimu1);
[results1, chain1, s2chain1] = mcmcrun(model, data, params, options);

save(fullfile(mcmc_dir, "Chain1_intermediate.mat"), 'results1', 'chain1', 's2chain1');
fprintf('Chain 1 checkpoint saved.\n');

options.nsimu      = nsimu2;
options.burnintime = burn2;

fprintf('Starting Chain 2 (continuation) (%g iterations)...\n', nsimu2);
[results, chain2, s2chain2] = mcmcrun(model, data, params, options, results1);

save(fullfile(mcmc_dir, "Chain2_final.mat"), 'results', 'chain2', 's2chain2');
fprintf('Chain 2 saved.\n');

% Post-burn samples for plotting/inference
burn = min(burn2, size(chain2,1)-1);
chain_post = chain2(burn+1:end, :);

%% [5] ----------------------- Publication-style plots --------------------
figHist  = plot_mcmcplot_2x4(chain_post, results, param_names, 'hist');
figTrace = plot_mcmcplot_2x4(chain_post, results, param_names, 'trace');

chainstats(chain_post, results);

%% [6] -------------------- Model fit using mean(theta) -------------------
y0 = [0.07035; 0; 3.40861; 0; 0.12131; 0.11329; 0.001316];
ode_options = odeset('NonNegative', 1:7, 'RelTol',1e-6, 'AbsTol',1e-9);

theta_mean = mean(chain_post, 1).';  % column vector
[t, y] = ode15s(@(t, y) Cancer4sys(t, y, theta_mean), data.xdata, y0, ode_options);

figFit = figure('Color','w'); clf(figFit);
plot(t, y(:,7), '-', 'LineWidth', 1.6, 'DisplayName', '$\mathrm{PI3K}$ (Model)');
hold on;
scatter(data.xdata, data.ydata, 35, 'filled', 'DisplayName', 'Data');
hold off;

xlabel('Time (minutes)', 'Interpreter','latex');
ylabel('Concentration ($\mu$M)', 'Interpreter','latex');
title('Model fit to $\mathrm{PI3K}$ data (post-burn mean)', 'Interpreter','latex');
legend('Location','southeast', 'Interpreter','latex', 'Box','off');
grid on; box on;
set(gca,'FontSize',10,'LineWidth',1);

%% [6.4] ---- Credible band (median + 95% CI) (separate fig) --------------
rng(1,'twister');
nsamp = 300;
tgrid = linspace(min(data.xdata), max(data.xdata), 200).';

[medI, loI, hiI] = posterior_band_PI3K(chain_post, tgrid, y0, ode_options, nsamp);

figBand = figure('Color','w'); clf(figBand);
ax = gca; hold on;
c  = ax.ColorOrder(1,:);

patch([tgrid; flipud(tgrid)], [loI; flipud(hiI)], c, ...
      'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility','off');

plot(tgrid, medI, '-', 'LineWidth', 1.8, 'Color', c, 'DisplayName', 'Posterior median');
scatter(data.xdata, data.ydata, 30, 'filled', 'DisplayName', 'Data');

xlabel('Time (minutes)', 'Interpreter','latex');
ylabel('Concentration ($\mu$M)', 'Interpreter','latex');
title('Credible band for $\mathrm{PI3K}$ (95\% credible interval)', 'Interpreter','latex');
legend('Location','southeast', 'Interpreter','latex', 'Box','off');
grid on; box on;
set(gca,'FontSize',10,'LineWidth',1);
hold off;

%% [7] ---------------------------- Save figures --------------------------
save_figure(figHist,  fullfile(mcmc_dir, "MCMC_Histograms_2x4_mcmcplotStyle"));
save_figure(figTrace, fullfile(mcmc_dir, "MCMC_Traces_2x4_mcmcplotStyle"));
save_figure(figFit,   fullfile(mcmc_dir, "Simulated_vs_Observed_PI3K_mean"));
save_figure(figBand,  fullfile(mcmc_dir, "CredibleBand_PI3K_95CI"));

%% [8] ------------------------ Save workspace and exit -------------------
close all;
save(fullfile(mcmc_dir, "MCMC_Workspace_" + init_time + ".mat"), '-v7.3');

disp('MCMC run complete.');
diary off;

%=========================================================================%
% Local helper functions
%=========================================================================%

function fig = plot_mcmcplot_2x4(chain, results, param_names, mode)
% 2x4 tiled panels, each panel generated using mcmcplot (code-1 style),
% then copied into a tiledlayout (code-2 organization).

    nPar  = size(chain,2);
    nRows = 2;
    nCols = 4;

    fig = figure('Color','w');
    set(fig,'Units','inches','Position',[1 1 12.0 6.0]);
    tl = tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact');

    for j = 1:nPar
        axDst = nexttile(tl);

        figTmp = figure('Visible','off', 'Color','w');

        switch lower(mode)
            case 'hist'
                mcmcplot(chain, j, results, 'hist');
            case 'trace'
                mcmcplot(chain, j, results);
            otherwise
                close(figTmp);
                error('Unknown mode. Use "hist" or "trace".');
        end

        axTmpAll = findobj(figTmp,'Type','axes');
        if isempty(axTmpAll)
            close(figTmp);
            continue;
        end

        nchild = arrayfun(@(a) numel(allchild(a)), axTmpAll);
        [~,k] = max(nchild);
        axTmp = axTmpAll(k);

        copyobj(allchild(axTmp), axDst);

        axDst.XLim   = axTmp.XLim;
        axDst.YLim   = axTmp.YLim;
        axDst.XScale = axTmp.XScale;
        axDst.YScale = axTmp.YScale;

        axDst.FontSize  = 9;
        axDst.LineWidth = 1;
        grid(axDst,'on'); box(axDst,'on');

        title(axDst, param2latex(param_names{j}), 'Interpreter','latex', 'FontSize',11);

        if strcmpi(mode,'trace')
            xlabel(axDst, 'Iteration', 'Interpreter','latex');
        end

        close(figTmp);
    end

    if strcmpi(mode,'hist')
        title(tl, 'Posterior histograms (mcmcplot style)', 'Interpreter','latex', 'FontSize',14);
    else
        title(tl, 'MCMC chain traces (mcmcplot style)', 'Interpreter','latex', 'FontSize',14);
    end
end

function s = param2latex(name)
% alpha_A -> \alpha_{A}, beta_3 -> \beta_{3}, K1 -> K_{1}

    name = string(name);

    name = replace(name, "alpha", "\alpha");
    name = replace(name, "beta",  "\beta");
    name = replace(name, "gamma", "\gamma");
    name = replace(name, "delta", "\delta");

    if contains(name, "_")
        parts = split(name, "_");
        base  = parts(1);
        sub   = join(parts(2:end), "_");
        s = "$" + base + "_{" + sub + "}$";
        return;
    end

    tok = regexp(name, "^(.*?)(\d+)$", "tokens", "once");
    if ~isempty(tok)
        s = "$" + tok{1} + "_{" + tok{2} + "}$";
    else
        s = "$" + name + "$";
    end
end

function save_figure(h, basepath)
% Robust save: accepts figure OR axes handle and saves .fig + vector .svg

    if isempty(h) || ~isgraphics(h)
        warning('save_figure:InvalidHandle','Invalid/empty handle, skipping save.');
        return;
    end

    if ~isgraphics(h, 'figure')
        hfig = ancestor(h, 'figure');
    else
        hfig = h;
    end

    if isempty(hfig) || ~isgraphics(hfig, 'figure')
        warning('save_figure:NoFigure','Could not resolve a figure handle, skipping save.');
        return;
    end

    basepath = string(basepath);
    figfile  = basepath + ".fig";
    svgfile  = basepath + ".svg";

    drawnow;

    savefig(hfig, char(figfile));

    try
        exportgraphics(hfig, char(svgfile), 'ContentType','vector');
    catch
        saveas(hfig, char(svgfile));
    end
end

function [medI, loI, hiI] = posterior_band_PI3K(chain_use, tgrid, y0, ode_options, nsamp)
% Draw nsamp posterior samples, simulate PI3K (state 7), return median and 95% band

    n = size(chain_use,1);
    nsamp = min(nsamp, n);

    idx = randi(n, nsamp, 1);
    Y = NaN(length(tgrid), nsamp);

    for k = 1:nsamp
        theta_k = chain_use(idx(k), :).';
        try
            [~, yk] = ode15s(@(t,y) Cancer4sys(t, y, theta_k), tgrid, y0, ode_options);
            Y(:,k) = yk(:,7);
        catch
            % keep NaN for failed samples
        end
    end

    good = all(isfinite(Y), 1);
    Y = Y(:,good);

    if isempty(Y)
        error('All posterior simulations failed. Reduce nsamp or check ODE stability.');
    end

    medI = median(Y, 2);
    loI  = prctile(Y.',  2.5).';
    hiI  = prctile(Y.', 97.5).';
end

