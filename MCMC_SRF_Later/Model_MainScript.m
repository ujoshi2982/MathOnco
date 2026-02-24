% Main script to solve the ODEs
%% ====== Data for x and y as column vectors (Time in Minutes and µM for concentration) ======
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
          
%%===================Define parameters=======================
% Define parameters
alpha_A = 1.3176;     % Basal rate of LYN synthesis (µM min^-1)
alpha_B = 0.0247;  % Basal rate of FcERI synthesis (µM min^-1)
alpha_C = 0.8967;     % Basal rate of FYN synthesis (µM min^-1)
alpha_D = 0.0550;     % Basal rate of SHIP synthesis (µM min^-1)

gamma_A = 2.1e-3;      % Degradation rate of LYN (min^-1)
gamma_B = 2.1e-3;   % Degradation rate of FcERI (min^-1)
gamma_C = 1.1e-3;     % Degradation rate of FYN (min^-1)
gamma_D = 1.1981;     % Degradation rate of Shp1 (min^-1)      
gamma_I = 1e-3;      % Degradation rate of PI3K (min^-1)

K1 = 1.77;       % Km for LYN_FcERI through SHIP (µM)
K2 = 1.4087;         % Km for FYN through LYN_FcERI (µM)
K3 = 2.3072;      % Km for PI3K through FYN_FcERI (µM)

beta_3 = 5.221;         % Vmax for activation of PI3K (µM min^-1)

Kf1 = 4.8;         % Association constant of LYN and FcERI (µM^-1 min^-1)
Kr1 = 0.6;         % Dissociation constant of LYN and FcERI (min^-1)
Kf2 = 4.8;         % Association constant of FYN and FcERI (µM^-1 min^-1)
Kr2 = 1.7438;         % Dissociation constant of FYN and FcERI (min^-1)


%Initial concentrations (µM)
A0 = 0.070366;          % Initial concentration of LYN
AB0 = 0;           % Initial concentration of LYN-FcERI complex
B0 = 3.40861;        % Initial concentration of FcERI
BC0 = 0;           % Initial concentration of FcERI-FYN complex
C0 = 0.121311;        % Initial concentration of FYN
D0 = 0.113287;        % Initial concentration of SHIP
I0 = 0.001316;         % Initial concentration of PI3K

%  A0 = 4.65e-2;          % Initial concentration of LYN
%  AB0 = 0;           % Initial concentration of LYN-FcERI complex
%  B0 = 0.664;        % Initial concentration of FcERI
%  BC0 = 0;           % Initial concentration of FcERI-FYN complex
%  C0 = 0.169;        % Initial concentration of FYN
%  D0 = 1.165;        % Initial concentration of SHIP
%  I0 = 0.01;         % Initial concentration of PI3K

%%

% Time span for the simulation
tspan = [0 30];   % Time in minutes

% Initial conditions
y0 = [A0, AB0, B0, BC0, C0, D0, I0];

% Solve the ODE system, pass parameters as additional arguments
[t, y] = ode15s(@(t, y) model(t, y, alpha_A, alpha_B, alpha_C, alpha_D, gamma_A, gamma_B, gamma_C, gamma_D, gamma_I, ...
                       Kf1, Kr1, Kf2, Kr2, K1, K2, K3, beta_3), tspan, y0);

% Plot results with bold lines of size 2.0
figure;
plot(t, y(:,1), 'b', 'LineWidth', 2.0, 'DisplayName', 'LYN'); hold on;
plot(t, y(:,2), 'r', 'LineWidth', 2.0, 'DisplayName', 'LYN-FcERI');
plot(t, y(:,3), 'g', 'LineWidth', 2.0, 'DisplayName', 'FcERI');
plot(t, y(:,4), 'c', 'LineWidth', 2.0, 'DisplayName', 'FcERI-FYN');
plot(t, y(:,5), 'm', 'LineWidth', 2.0, 'DisplayName', 'FYN');
plot(t, y(:,6), 'Color', [0.929 0.694 0.125], 'LineWidth', 2.0, 'DisplayName', 'SHIP');
plot(t, y(:,7), 'Color', [0.494 0.184 0.556], 'LineWidth', 2.0, 'DisplayName', 'PI3K');
scatter(data.xdata, data.ydata, 50, 'r', 'DisplayName', 'Observed Data');
hold off;

% Set labels and title
xlabel('Time (minutes)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Concentration (µM)', 'FontWeight', 'bold', 'FontSize', 12);
% Set limits for x-axis and y-axis
%xlim([0 1e3]); % Example: X-axis from 0 to 100
%ylim([0 30]);  % Example: Y-axis from 0 to 50
title('SHIP, Con=0, basal, deg=zero', 'FontWeight', 'bold', 'FontSize', 14);

% Display legend
legend('show', 'Location', 'northeast');

% Add grid
grid on;