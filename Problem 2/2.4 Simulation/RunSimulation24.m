
%% ---------------------------------------------------------
% Problem 2.4 - Open-Loop Simulation
% Modified Four-Tank System
%
% Models:
% 2.1 Deterministic nonlinear model
% 2.2 Piecewise-constant stochastic disturbances
% 2.3 Ornstein-Uhlenbeck SDE model
%
% Dynamics:     dx/dt = f(t,x,u,d,p)
% Measurements: y = g(x,p) + v
% Outputs:      z = h(x,p)
% ---------------------------------------------------------

clear; clc; close all;

%% ---------------------------------------------------------
% Add function folders to MATLAB path
% ---------------------------------------------------------

scriptFolder = fileparts(mfilename('fullpath'));
parentFolder = fileparts(scriptFolder);

addpath(fullfile(parentFolder,'Common'));
addpath(fullfile(parentFolder,'2.3 SDE Model'));

%% ---------------------------------------------------------
% Physical parameters
% ---------------------------------------------------------

a = 1.2272*ones(4,1);
A = 380.1327*ones(4,1);

gamma = [0.6;0.7];

g = 981;
rho = 1.00;

p = [a;A;gamma;g;rho];

%% ---------------------------------------------------------
% Simulation settings
% ---------------------------------------------------------

tf = 20*60;
Ts = 4;
dt = 0.25;

T = (0:Ts:tf)';
N = length(T);

%% ---------------------------------------------------------
% Piecewise-constant manipulated inputs
% ---------------------------------------------------------

% Each row contains [F1 F2] in cm^3/s
Usteps = [300 300;
          330 300;
          330 330;
          300 330];

% Switching times [s]
switchTimes = [0 300 600 900 1200];

% Generate manipulated input sequence
U = zeros(N,2);

for k = 1:N

    tk = T(k);

    if tk < 300
        U(k,:) = Usteps(1,:);
    elseif tk < 600
        U(k,:) = Usteps(2,:);
    elseif tk < 900
        U(k,:) = Usteps(3,:);
    else
        U(k,:) = Usteps(4,:);
    end

end

%% ---------------------------------------------------------
% Disturbance parameters
% ---------------------------------------------------------

% Mean disturbances for all models
dbar = [20;20];

% Piecewise-constant stochastic model
d_std = [5;5];
Td = 20;

% Ornstein-Uhlenbeck model
ad = [1/60;1/60];
stdF = [5;5];

sigmaF = stdF.*sqrt(2*ad);

%% ---------------------------------------------------------
% Measurement noise parameters
% ---------------------------------------------------------

sigma_y = [0.1;0.1];

Rvv = diag(sigma_y.^2);

Lv = chol(Rvv,'lower');

%% ---------------------------------------------------------
% Common initial conditions
% Nominal steady state
% ---------------------------------------------------------

u0 = [300;300];

q3s = (1-gamma(2))*u0(2) + dbar(1);
q4s = (1-gamma(1))*u0(1) + dbar(2);

q1s = gamma(1)*u0(1) + q3s;
q2s = gamma(2)*u0(2) + q4s;

qs = [q1s;q2s;q3s;q4s];

hs = (qs./a).^2/(2*g);

m0 = rho*A.*hs;

%% =========================================================
% MODEL 2.1 - Deterministic nonlinear system
% ==========================================================

Xdet = zeros(N,4);
Xdet(1,:) = m0';

for k = 1:N-1

    uk = U(k,:)';
    dk = dbar;

    [~,Xseg] = ode45( ...
        @(tt,xx) FourTankSystem(tt,xx,uk,dk,p), ...
        [T(k) T(k+1)],Xdet(k,:)');

    Xdet(k+1,:) = Xseg(end,:);

end

%% =========================================================
% MODEL 2.2 - Piecewise-constant stochastic disturbances
% ==========================================================

rng(1);

Xpc = zeros(N,4);
Xpc(1,:) = m0';

Dpc = zeros(N,2);

assert(mod(tf,Td)==0, ...
    'tf must be a multiple of Td');

assert(mod(Td,Ts)==0, ...
    'Td must be a multiple of Ts');

nDist = tf/Td;

Dvalues = dbar + d_std.*randn(2,nDist);

% Keep inlet disturbances nonnegative
Dvalues = max(Dvalues,0);

for k = 1:N-1

    uk = U(k,:)';

    interval = floor(T(k)/Td)+1;
    dk = Dvalues(:,interval);

    Dpc(k,:) = dk';

    [~,Xseg] = ode45( ...
        @(tt,xx) FourTankSystem(tt,xx,uk,dk,p), ...
        [T(k) T(k+1)],Xpc(k,:)');

    Xpc(k+1,:) = Xseg(end,:);

end

Dpc(end,:) = Dpc(end-1,:);

%% =========================================================
% MODEL 2.3 - Ornstein-Uhlenbeck SDE
% ==========================================================

Nfine = round(tf/dt);

Tf = (0:Nfine)*dt;

% Two independent Wiener processes
rng(100);
dW = sqrt(dt)*randn(2,Nfine);

% Augmented state [m1;m2;m3;m4;F3;F4]
Xsde = zeros(6,Nfine+1);

Xsde(:,1) = [m0;dbar];

for k = 1:Nfine

    tk = Tf(k);

    % Piecewise-constant pump inputs
    if tk < 300
        uk = Usteps(1,:)';
    elseif tk < 600
        uk = Usteps(2,:)';
    elseif tk < 900
        uk = Usteps(3,:)';
    else
        uk = Usteps(4,:)';
    end

    xk = Xsde(:,k);

    % Drift
    f = FourTankDrift( ...
        tk,xk,uk,p,ad,dbar,sigmaF);

    % Diffusion
    sigma = FourTankDiffusion( ...
        tk,xk,uk,p,ad,dbar,sigmaF);

    % Euler-Maruyama step
    Xsde(:,k+1) = xk + f*dt + sigma*dW(:,k);

end

%% ---------------------------------------------------------
% Sample SDE at measurement time Ts
% ---------------------------------------------------------

assert(abs(Ts/dt-round(Ts/dt))<1e-10, ...
    'Ts must be a multiple of dt');

sampleIndex = 1:round(Ts/dt):Nfine+1;

XsdeSample = Xsde(:,sampleIndex);

Dsde = XsdeSample(5:6,:)';

%% ---------------------------------------------------------
% Check numerical feasibility
% ---------------------------------------------------------

if any(~isfinite(Xdet(:))) || ...
   any(~isfinite(Xpc(:))) || ...
   any(~isfinite(Xsde(:)))

    error('Simulation generated non-finite states.');

end

if any(Xdet(:)<0) || any(Xpc(:)<0) || ...
   any(any(Xsde(1:4,:)<0))

    error('Negative tank masses encountered.');

end

if any(any(Xsde(5:6,:)<0))
    warning('Negative OU disturbance flows encountered.');
end

%% ---------------------------------------------------------
% Controlled output and measurement models
%
% z(t_k) = h(x(t_k),p)
% y(t_k) = g(x(t_k),p) + v_k
% ---------------------------------------------------------

Zdet = zeros(N,2);
Zpc  = zeros(N,2);
Zsde = zeros(N,2);

Ydet = zeros(N,2);
Ypc  = zeros(N,2);
Ysde = zeros(N,2);

% Reproducible sensor noise
rng(200);

Vdet = randn(N,2)*Lv';
Vpc  = randn(N,2)*Lv';
Vsde = randn(N,2)*Lv';

for k = 1:N

    % Physical state vectors
    xdet = Xdet(k,:)';
    xpc  = Xpc(k,:)';
    xsde = XsdeSample(1:4,k);

    % Controlled outputs
    Zdet(k,:) = FourTankOutput(xdet,p)';
    Zpc(k,:)  = FourTankOutput(xpc,p)';
    Zsde(k,:) = FourTankOutput(xsde,p)';

    % Sensor measurements with Gaussian noise
    Ydet(k,:) = FourTankSensor(xdet,p)' + Vdet(k,:);
    Ypc(k,:)  = FourTankSensor(xpc,p)'  + Vpc(k,:);
    Ysde(k,:) = FourTankSensor(xsde,p)' + Vsde(k,:);

end

%% =========================================================
% FIGURE 1 - Manipulated inputs
% ==========================================================

figure;

subplot(2,1,1)

stairs(T/60,U(:,1),'LineWidth',1.5)
xlabel('Time [min]')
ylabel('F_1 [cm^3/s]')
title('Manipulated Input F1')
grid on

subplot(2,1,2)

stairs(T/60,U(:,2),'LineWidth',1.5)
xlabel('Time [min]')
ylabel('F_2 [cm^3/s]')
title('Manipulated Input F2')
grid on

sgtitle('Open-Loop Piecewise-Constant Inputs')

%% =========================================================
% FIGURE 2 - Disturbance trajectories
% ==========================================================

figure;

for j = 1:2

    subplot(2,1,j)
    hold on

    stairs(T/60,Dpc(:,j), ...
        'LineWidth',1.2, ...
        'DisplayName','Piecewise constant');

    plot(T/60,Dsde(:,j), ...
        'LineWidth',1.2, ...
        'DisplayName','OU process');

    yline(dbar(j),'k--', ...
        'LineWidth',1.2, ...
        'DisplayName','Mean disturbance');

    xlabel('Time [min]')
    ylabel(sprintf('F_%d [cm^3/s]',j+2))
    title(sprintf('Disturbance F%d',j+2))

    legend('Location','best')
    grid on
    hold off

end

sgtitle('Unmeasured Disturbance Trajectories')

%% =========================================================
% FIGURE 3 - Controlled output comparison
% ==========================================================

figure;

for j = 1:2

    subplot(2,1,j)
    hold on

    plot(T/60,Zdet(:,j),'k', ...
        'LineWidth',1.8, ...
        'DisplayName','Deterministic');

    plot(T/60,Zpc(:,j),'--', ...
        'LineWidth',1.5, ...
        'DisplayName','Piecewise stochastic');

    plot(T/60,Zsde(:,j),':', ...
        'LineWidth',1.8, ...
        'DisplayName','OU SDE');

    xlabel('Time [min]')
    ylabel(sprintf('h_%d [cm]',j))
    title(sprintf('Tank %d - Controlled Output',j))

    legend('Location','best')
    grid on
    hold off

end

sgtitle('Open-Loop Output Comparison')

%% =========================================================
% FIGURE 4 - True outputs vs noisy measurements
% OU SDE model
% ==========================================================

figure;

for j = 1:2

    subplot(2,1,j)
    hold on

    plot(T/60,Zsde(:,j),'b', ...
        'LineWidth',1.5, ...
        'DisplayName','True output');

    plot(T/60,Ysde(:,j),'r.', ...
        'MarkerSize',5, ...
        'DisplayName','Noisy measurement');

    xlabel('Time [min]')
    ylabel(sprintf('h_%d [cm]',j))
    title(sprintf('Tank %d',j))

    legend('Location','best')
    grid on
    hold off

end

sgtitle('OU Model - Outputs and Measurements')

%% ---------------------------------------------------------
% Display simulation information
% ---------------------------------------------------------

disp('-----------------------------------------')
disp('Problem 2.4 - Open-Loop Simulation')
disp('-----------------------------------------')

disp('Manipulated input sequence [F1 F2]:')
disp(Usteps)

disp('Initial steady-state levels [cm]:')
disp(hs)

disp('Measurement covariance Rvv [cm^2]:')
disp(Rvv)

fprintf('Simulation time = %.0f s\n',tf);
fprintf('Sampling time = %.1f s\n',Ts);
fprintf('SDE integration step = %.2f s\n',dt);
fprintf('Disturbance update interval = %.0f s\n',Td);
