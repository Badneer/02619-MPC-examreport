
%% ---------------------------------------------------------
% Problem 2.3 - Stochastic Nonlinear Model
% Modified Four-Tank System
% Ornstein-Uhlenbeck disturbance model
% Euler-Maruyama numerical method
% ---------------------------------------------------------

clear; clc; close all;

%% Add Common folder to path

scriptFolder = fileparts(mfilename('fullpath'));
problem2Folder = fileparts(scriptFolder);

addpath(fullfile(problem2Folder,'Common'));

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
% Manipulated variables
% ---------------------------------------------------------

F1 = 300;
F2 = 300;

u = [F1;F2];

%% ---------------------------------------------------------
% Ornstein-Uhlenbeck disturbance parameters
% ---------------------------------------------------------

% Mean disturbance flows [cm^3/s]
dbar = [20;20];

% Mean-reversion rates [1/s]
ad = [1/60;1/60];

% Desired stationary standard deviations [cm^3/s]
stdF = [5;5];

% Diffusion coefficients
sigmaF = stdF.*sqrt(2*ad);

%% ---------------------------------------------------------
% Simulation settings
% ---------------------------------------------------------

tf = 20*60;
dt = 0.25;
Ts = 4;

N = round(tf/dt);

nW = 2;
Ns = 5;
seed = 100;

assert(abs(tf-N*dt)<1e-10, ...
    'tf must be a multiple of dt');

assert(abs(Ts/dt-round(Ts/dt))<1e-10, ...
    'Ts must be a multiple of dt');

%% ---------------------------------------------------------
% Initial conditions
% Nominal steady state
% ---------------------------------------------------------

q3s = (1-gamma(2))*F2 + dbar(1);
q4s = (1-gamma(1))*F1 + dbar(2);

q1s = gamma(1)*F1 + q3s;
q2s = gamma(2)*F2 + q4s;

qs = [q1s;q2s;q3s;q4s];

hs = (qs./a).^2/(2*g);

m0 = rho*A.*hs;

% Augmented initial state: 4 masses + 2 flows
x0 = [m0;dbar];

%% ---------------------------------------------------------
% Generate Wiener processes
% ---------------------------------------------------------

[W,T,~] = StdWienerProcess( ...
    tf,N,nW,Ns,seed);

%% ---------------------------------------------------------
% SDE simulation
% ---------------------------------------------------------

X = zeros(6,N+1,Ns);

for i = 1:Ns

    X(:,:,i) = SDEsolverExplicitExplicit( ...
        @FourTankDrift, ...
        @FourTankDiffusion, ...
        T,x0,W(:,:,i), ...
        u,p,ad,dbar,sigmaF);

end

%% ---------------------------------------------------------
% Check numerical and physical feasibility
% ---------------------------------------------------------

if any(~isfinite(X(:)))
    error('Non-finite states encountered.');
end

masses = X(1:4,:,:);
disturbances = X(5:6,:,:);

if any(masses(:)<0)
    error(['Negative masses generated. ', ...
        'The simulation/model requires review.']);
end

if any(disturbances(:)<0)
    warning(['Negative disturbance flow encountered. ', ...
        'OU processes do not enforce positivity.']);
end

%% ---------------------------------------------------------
% Calculate physical tank levels
% ---------------------------------------------------------

H = zeros(4,N+1,Ns);

for i = 1:Ns

    H(:,:,i) = X(1:4,:,i)./(rho*A);

end

%% ---------------------------------------------------------
% Measurement and output models
% ---------------------------------------------------------

% Measurement covariance
Rvv = diag([0.1^2;0.1^2]);

Lv = chol(Rvv,'lower');

% Measurement sampling points
sampleIndex = 1:round(Ts/dt):N+1;

Ty = T(sampleIndex);
Ny = length(sampleIndex);

Y = zeros(2,Ny,Ns);
Z = zeros(2,Ny,Ns);

for i = 1:Ns

    for k = 1:Ny

        idx = sampleIndex(k);

        % Only four physical mass states
        xphysical = X(1:4,idx,i);

% Calculate true liquid levels of Tanks 1 and 2
zk = xphysical(1:2)./(rho*A(1:2));

% True controlled outputs
Z(:,k,i) = zk;

% Gaussian measurement noise
vk = Lv*randn(2,1);

% Noisy sensor measurements
Y(:,k,i) = zk + vk;


    end

end

%% ---------------------------------------------------------

%% ---------------------------------------------------------
% Figure 1 - Stochastic disturbance flows
% ---------------------------------------------------------

figure;

realization_labels = arrayfun( ...
    @(i) sprintf('Realization %d',i), ...
    1:Ns,'UniformOutput',false);

for j = 1:2

    subplot(2,1,j)
    hold on

    for i = 1:Ns

        plot(T/60, ...
            squeeze(X(j+4,:,i)), ...
            'LineWidth',1.2, ...
            'DisplayName',realization_labels{i});

    end

    % Mean disturbance flow
    yline(dbar(j),'k--', ...
        'LineWidth',1.5, ...
        'DisplayName','Mean flow');

    xlabel('Time [min]')
    ylabel(sprintf('F_%d [cm^3/s]',j+2))
    title(sprintf('Disturbance F%d',j+2))

    legend('Location','best','FontSize',8)

    grid on
    hold off

end

sgtitle('Ornstein-Uhlenbeck Disturbances')



%% ---------------------------------------------------------
% Figure 2 - Tank masses
% ---------------------------------------------------------

figure;

tank_order = [3 4 1 2];

for j = 1:4

    tank = tank_order(j);

    subplot(2,2,j)
    hold on

    for i = 1:Ns

        plot(T/60, ...
            squeeze(X(tank,:,i))/1000, ...
            'LineWidth',1.2, ...
            'DisplayName',realization_labels{i});

    end

    xlabel('Time [min]')
    ylabel(sprintf('m_%d [kg]',tank))
    title(sprintf('Tank %d',tank))

    legend('Location','best','FontSize',8)

    grid on
    hold off

end

sgtitle('Stochastic Tank Masses')


%% ---------------------------------------------------------
% Figure 3 - Tank levels
% ---------------------------------------------------------

figure;

for j = 1:4

    tank = tank_order(j);

    subplot(2,2,j)
    hold on

    for i = 1:Ns

        plot(T/60, ...
            squeeze(H(tank,:,i)), ...
            'LineWidth',1.2, ...
            'DisplayName',realization_labels{i});

    end

    % Nominal steady-state level
    yline(hs(tank),'k--', ...
        'LineWidth',1.5, ...
        'DisplayName','Nominal steady state');

    xlabel('Time [min]')
    ylabel(sprintf('h_%d [cm]',tank))
    title(sprintf('Tank %d',tank))

    legend('Location','best','FontSize',8)

    grid on
    hold off

end

sgtitle('Stochastic Tank Levels')


%% ---------------------------------------------------------
% Figure 4 - Measurements
% First realization
% ---------------------------------------------------------

figure;

for j = 1:2

    subplot(2,1,j)

    plot(Ty/60,squeeze(Z(j,:,1)), ...
        'b','LineWidth',1.5)

    hold on

    plot(Ty/60,squeeze(Y(j,:,1)), ...
        'r.','MarkerSize',5)

    xlabel('Time [min]')
    ylabel(sprintf('h_%d [cm]',j))
    title(sprintf('Tank %d',j))

    legend('True output','Noisy measurement', ...
        'Location','best')

    grid on

end

sgtitle('Controlled Outputs and Measurements')

%% ---------------------------------------------------------
% Display results
% ---------------------------------------------------------

disp('-----------------------------------------')
disp('Problem 2.3 - OU Stochastic Model')
disp('-----------------------------------------')

disp('Nominal pump flows [cm^3/s]:')
disp(u)

disp('Mean disturbance flows [cm^3/s]:')
disp(dbar)

disp('Mean-reversion coefficients [1/s]:')
disp(ad)

disp('Diffusion coefficients:')
disp(sigmaF)

disp('Steady-state levels [cm]:')
disp(hs)

disp('Measurement covariance Rvv [cm^2]:')
disp(Rvv)

fprintf('Integration time step = %.3f s\n',dt);
fprintf('Measurement sampling time = %.1f s\n',Ts);
fprintf('Number of realizations = %d\n',Ns);
