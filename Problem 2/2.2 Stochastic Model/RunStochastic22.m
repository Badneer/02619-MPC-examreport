
%% ---------------------------------------------------------
% Problem 2.2 - Stochastic Nonlinear Model
% Modified Four-Tank System
%
% Dynamics:     dx/dt = f(t,x,u,d,p)
% Measurements: y = g(x,p) + v
% Outputs:      z = h(x,p)
% ---------------------------------------------------------

clear; clc; close all;

rng(1);  % Reproducible random numbers

%% ---------------------------------------------------------
% Add Common folder to MATLAB path
% ---------------------------------------------------------

scriptFolder = fileparts(mfilename('fullpath'));
parentFolder = fileparts(scriptFolder);

addpath(fullfile(parentFolder,'Common'));

%% ---------------------------------------------------------
% Physical parameters (same as Problem 2.1)
% ---------------------------------------------------------

a1 = 1.2272;
a2 = 1.2272;
a3 = 1.2272;
a4 = 1.2272;

A1 = 380.1327;
A2 = 380.1327;
A3 = 380.1327;
A4 = 380.1327;

gamma1 = 0.6;
gamma2 = 0.7;

g = 981;
rho = 1.00;

p = [a1;a2;a3;a4; A1;A2;A3;A4; ...
     gamma1;gamma2;g;rho];

A = [A1;A2;A3;A4];

%% ---------------------------------------------------------
% Simulation scenario
% ---------------------------------------------------------

t0 = 0;
tf = 20*60;

Ts = 4;    % Measurement sampling time [s]
Td = 20;   % Disturbance update interval [s]

t = (t0:Ts:tf)';
N = length(t);

% Initial conditions: empty tanks
x0 = [0;0;0;0];

%% ---------------------------------------------------------
% Nominal manipulated variables
% ---------------------------------------------------------

F1_nom = 300;   % Pump 1 [cm^3/s]
F2_nom = 300;   % Pump 2 [cm^3/s]

u = [F1_nom;F2_nom];

%% ---------------------------------------------------------
% Stochastic disturbance parameters
% ---------------------------------------------------------

% Mean disturbance flows [cm^3/s]
d_mean = [20;20];

% Standard deviations [cm^3/s]
d_std = [5;5];

% Disturbance covariance matrix
Qd = diag(d_std.^2);

%% ---------------------------------------------------------
% Measurement noise parameters
% ---------------------------------------------------------

% Sensor standard deviations [cm]
sigma_y = [0.1;0.1];

% Measurement noise covariance matrix
Rvv = diag(sigma_y.^2);

%% ---------------------------------------------------------
% Storage
% ---------------------------------------------------------

X = zeros(N,4);
D = zeros(N,2);

% Initial condition
X(1,:) = x0';

%% ---------------------------------------------------------
% Stochastic nonlinear simulation
% ---------------------------------------------------------

% Number of sampling steps per disturbance interval
stepsPerInterval = round(Td/Ts);

% Verify interval compatibility
assert(abs(stepsPerInterval*Ts-Td)<1e-10, ...
    'Td must be a multiple of Ts');

assert(mod(N-1,stepsPerInterval)==0, ...
    'Simulation time must contain complete intervals');

nIntervals = (N-1)/stepsPerInterval;

for k = 1:nIntervals

    % Random disturbance values
    dk = d_mean + d_std.*randn(2,1);

    % Prevent negative disturbance inlet flows
    dk = max(dk,0);

    % Time indices for current disturbance interval
    i0 = (k-1)*stepsPerInterval + 1;
    i1 = k*stepsPerInterval + 1;

    % Solve nonlinear dynamics with constant disturbance
    [~,Xseg] = ode15s( ...
        @(tt,xx) FourTankSystem(tt,xx,u,dk,p), ...
        t(i0:i1),X(i0,:)');

    % Store simulated tank masses
    X(i0:i1,:) = Xseg;

    % Store piecewise-constant disturbances
    D(i0:i1-1,:) = repmat( ...
        dk',stepsPerInterval,1);

end

D(end,:) = D(end-1,:);

%% ---------------------------------------------------------
% Calculate all four tank levels
% ---------------------------------------------------------

H = X./(rho*A');

%% ---------------------------------------------------------
% Controlled output and measurement models
%
% z(t_k) = h(x(t_k),p)
% y(t_k) = g(x(t_k),p) + v_k
% ---------------------------------------------------------

Z = zeros(N,2);
Y = zeros(N,2);

% Generate Gaussian measurement noise
Lv = chol(Rvv,'lower');

V = randn(N,2)*Lv';

for k = 1:N

    % Physical state vector [m1;m2;m3;m4]
    xk = X(k,:)';

    % Controlled outputs (true levels)
    Z(k,:) = FourTankOutput(xk,p)';

    % Noise-free sensor measurements
    gk = FourTankSensor(xk,p);

    % Add Gaussian measurement noise
    Y(k,:) = gk' + V(k,:);

end

%% ---------------------------------------------------------
% Plot 1 - Tank masses
% ---------------------------------------------------------

figure;

tank_order = [3 4 1 2];

for j = 1:4

    tank = tank_order(j);

    subplot(2,2,j)

    plot(t/60,X(:,tank)/1000, ...
        'b','LineWidth',1.5)

    xlabel('Time [min]')
    ylabel(sprintf('m_%d [kg]',tank))
    title(sprintf('Tank %d',tank))

    xlim([0 tf/60])
    grid on

end

sgtitle('Stochastic Nonlinear Model - Tank Masses')

%% ---------------------------------------------------------
% Plot 2 - Stochastic disturbances
% ---------------------------------------------------------

figure;

subplot(2,1,1)

stairs(t/60,D(:,1),'LineWidth',1.5)

xlabel('Time [min]')
ylabel('F_3 [cm^3/s]')
title('Stochastic Disturbance F3')
xlim([0 tf/60])
grid on

subplot(2,1,2)

stairs(t/60,D(:,2),'LineWidth',1.5)

xlabel('Time [min]')
ylabel('F_4 [cm^3/s]')
title('Stochastic Disturbance F4')
xlim([0 tf/60])
grid on

sgtitle('Piecewise-Constant Stochastic Disturbances')

%% ---------------------------------------------------------
% Plot 3 - True levels vs noisy measurements
% ---------------------------------------------------------

figure;

subplot(2,1,1)

plot(t/60,Z(:,1),'b','LineWidth',1.5)
hold on

plot(t/60,Y(:,1),'r.','MarkerSize',5)

xlabel('Time [min]')
ylabel('h_1 [cm]')
title('Tank 1 - True Level vs Measurement')

legend('True level','Noisy measurement', ...
    'Location','best')

xlim([0 tf/60])
grid on

subplot(2,1,2)

plot(t/60,Z(:,2),'b','LineWidth',1.5)
hold on

plot(t/60,Y(:,2),'r.','MarkerSize',5)

xlabel('Time [min]')
ylabel('h_2 [cm]')
title('Tank 2 - True Level vs Measurement')

legend('True level','Noisy measurement', ...
    'Location','best')

xlim([0 tf/60])
grid on

sgtitle('Controlled Outputs and Measurements')

%% ---------------------------------------------------------
% Plot 4 - All tank levels
% ---------------------------------------------------------

figure;

for j = 1:4

    tank = tank_order(j);

    subplot(2,2,j)

    plot(t/60,H(:,tank), ...
        'b','LineWidth',1.5)

    xlabel('Time [min]')
    ylabel(sprintf('h_%d [cm]',tank))
    title(sprintf('Tank %d',tank))

    xlim([0 tf/60])
    grid on

end

sgtitle('Stochastic Nonlinear Model - Tank Levels')

%% ---------------------------------------------------------
% Display simulation information
% ---------------------------------------------------------

disp('---------------------------------------')
disp('Problem 2.2 - Stochastic Nonlinear Model')
disp('---------------------------------------')

fprintf('Nominal F1 = %.1f cm^3/s\n',F1_nom)
fprintf('Nominal F2 = %.1f cm^3/s\n',F2_nom)

fprintf('Mean F3 = %.1f cm^3/s\n',d_mean(1))
fprintf('Mean F4 = %.1f cm^3/s\n',d_mean(2))

fprintf('Disturbance interval = %.1f s\n',Td)
fprintf('Sampling time = %.1f s\n',Ts)

disp('Disturbance covariance matrix Qd:')
disp(Qd)

disp('Measurement covariance matrix Rvv:')
disp(Rvv)
