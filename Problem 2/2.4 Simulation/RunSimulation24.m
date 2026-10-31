
%% ---------------------------------------------------------
% Problem 2.4 - Open-Loop Simulation
% Modified Four-Tank System
%
% Compare:
% 2.1 Deterministic nonlinear model
% 2.2 Piecewise-constant stochastic disturbances
% 2.3 Ornstein-Uhlenbeck SDE model
% ---------------------------------------------------------

clear; clc; close all;

%% ---------------------------------------------------------
% Add function folders to MATLAB path
% ---------------------------------------------------------

scriptFolder = fileparts(mfilename('fullpath'));
problem2Folder = fileparts(scriptFolder);

addpath(fullfile(problem2Folder,'Common'));
addpath(fullfile(problem2Folder,'2.3 SDE Model'));

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

tf = 20*60;       % Total time [s]
Ts = 4;           % Output sampling interval [s]
dt = 0.25;        % SDE integration step [s]

T = (0:Ts:tf)';
N = length(T);

%% ---------------------------------------------------------
% Piecewise-constant manipulated inputs
% ---------------------------------------------------------

% Rows: time intervals
% Columns: F1 and F2

Usteps = [300 300;
          330 300;
          330 330;
          300 330];

% Input switching times [s]
switchTimes = [0 300 600 900 1200];

% Generate inputs at the output sampling times
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

% Same mean disturbances for all three models
dbar = [20;20];       % [cm^3/s]

% Problem 2.2 parameters
d_std = [5;5];        % [cm^3/s]
Td = 20;              % Disturbance interval [s]

% Problem 2.3 OU parameters
ad = [1/60;1/60];     % Mean-reversion rates [1/s]
stdF = [5;5];         % Stationary std [cm^3/s]

sigmaF = stdF.*sqrt(2*ad);

%% ---------------------------------------------------------
% Common initial condition
% Steady state at F1 = F2 = 300 and d = dbar
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

    % Constant mean disturbances
    dk = dbar;

    [~,Xseg] = ode45( ...
        @(t,x) FourTankSystem(t,x,uk,dk,p), ...
        [T(k) T(k+1)],Xdet(k,:)');

    Xdet(k+1,:) = Xseg(end,:);

end

% True controlled outputs
Zdet = Xdet(:,1:2)./(rho*A(1:2)');

%% =========================================================
% MODEL 2.2 - Piecewise-constant stochastic model
% ==========================================================

rng(1);

Xpc = zeros(N,4);
Xpc(1,:) = m0';

Dpc = zeros(N,2);

% One disturbance value every Td seconds
nDist = tf/Td;

assert(mod(nDist,1)==0, ...
    'tf must be a multiple of Td');

Dvalues = dbar + d_std.*randn(2,nDist);

% Prevent negative disturbance inlet flows
Dvalues = max(Dvalues,0);

for k = 1:N-1

    uk = U(k,:)';

    % Find current disturbance interval
    interval = floor(T(k)/Td)+1;

    dk = Dvalues(:,interval);

    Dpc(k,:) = dk';

    % Integrate over one sampling interval
    [~,Xseg] = ode45( ...
        @(t,x) FourTankSystem(t,x,uk,dk,p), ...
        [T(k) T(k+1)],Xpc(k,:)');

    Xpc(k+1,:) = Xseg(end,:);

end

Dpc(end,:) = Dpc(end-1,:);

% True controlled outputs
Zpc = Xpc(:,1:2)./(rho*A(1:2)');

%% =========================================================
% MODEL 2.3 - Ornstein-Uhlenbeck SDE
% ==========================================================

Nfine = round(tf/dt);

% Fine integration grid
Tf = (0:Nfine)*dt;

% Wiener increments for two disturbance processes
rng(100);
dW = sqrt(dt)*randn(2,Nfine);

% Augmented state:
% [m1;m2;m3;m4;F3;F4]

Xsde = zeros(6,Nfine+1);
Xsde(:,1) = [m0;dbar];

for k = 1:Nfine

    tk = Tf(k);

    % Piecewise-constant manipulated input
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

    % Drift function
    f = FourTankDrift( ...
        tk,xk,uk,p,ad,dbar,sigmaF);

    % Diffusion function
    sigma = FourTankDiffusion( ...
        tk,xk,uk,p,ad,dbar,sigmaF);

    % Euler-Maruyama integration
    Xsde(:,k+1) = xk + f*dt + sigma*dW(:,k);

end

%% ---------------------------------------------------------
% Sample the SDE at Ts = 4 seconds
% ---------------------------------------------------------

sampleIndex = 1:round(Ts/dt):Nfine+1;

XsdeSample = Xsde(:,sampleIndex);

% True outputs
Zsde = XsdeSample(1:2,:)'./(rho*A(1:2)');

% OU disturbance flows
Dsde = XsdeSample(5:6,:)';

%% ---------------------------------------------------------
% Check numerical and physical feasibility
% ---------------------------------------------------------

allMasses = [Xdet(:); Xpc(:)];

if any(~isfinite(allMasses)) || any(~isfinite(Xsde(:)))
    error('Simulation generated non-finite states.');
end

if any(Xdet(:)<0) || any(Xpc(:)<0) || ...
   any(any(Xsde(1:4,:)<0))

    error('Negative tank masses encountered.');

end

if any(any(Xsde(5:6,:)<0))
    warning('OU disturbance flows became negative.');
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
% FIGURE 2 - Disturbances
% ==========================================================

figure;

subplot(2,1,1)

stairs(T/60,Dpc(:,1),'LineWidth',1.2)
hold on
plot(T/60,Dsde(:,1),'LineWidth',1.2)
yline(dbar(1),'k--','LineWidth',1.2)

xlabel('Time [min]')
ylabel('F_3 [cm^3/s]')
title('Disturbance F3')
legend('Piecewise constant','OU process', ...
       'Mean','Location','best')
grid on

subplot(2,1,2)

stairs(T/60,Dpc(:,2),'LineWidth',1.2)
hold on
plot(T/60,Dsde(:,2),'LineWidth',1.2)
yline(dbar(2),'k--','LineWidth',1.2)

xlabel('Time [min]')
ylabel('F_4 [cm^3/s]')
title('Disturbance F4')
legend('Piecewise constant','OU process', ...
       'Mean','Location','best')
grid on

sgtitle('Unmeasured Disturbance Trajectories')

%% =========================================================
% FIGURE 3 - Compare controlled outputs
% ==========================================================

figure;

subplot(2,1,1)

plot(T/60,Zdet(:,1),'k','LineWidth',1.8)
hold on
plot(T/60,Zpc(:,1),'--','LineWidth',1.5)
plot(T/60,Zsde(:,1),':','LineWidth',1.8)

xlabel('Time [min]')
ylabel('h_1 [cm]')
title('Tank 1 - Controlled Output')
legend('Deterministic','Piecewise stochastic', ...
       'OU SDE','Location','best')
grid on

subplot(2,1,2)

plot(T/60,Zdet(:,2),'k','LineWidth',1.8)
hold on
plot(T/60,Zpc(:,2),'--','LineWidth',1.5)
plot(T/60,Zsde(:,2),':','LineWidth',1.8)

xlabel('Time [min]')
ylabel('h_2 [cm]')
title('Tank 2 - Controlled Output')
legend('Deterministic','Piecewise stochastic', ...
       'OU SDE','Location','best')
grid on

sgtitle('Open-Loop Output Comparison')

%% ---------------------------------------------------------
% Print simulation settings
% ---------------------------------------------------------

disp('-----------------------------------------')
disp('Problem 2.4 - Open-Loop Simulation')
disp('-----------------------------------------')

disp('Pump input sequence [F1 F2]:')
disp(Usteps)

disp('Initial nominal levels [cm]:')
disp(hs)

fprintf('Simulation time: %.0f s\n',tf);
fprintf('Output sampling interval: %.1f s\n',Ts);
fprintf('SDE integration step: %.2f s\n',dt);
fprintf('Piecewise disturbance interval: %.0f s\n',Td);
