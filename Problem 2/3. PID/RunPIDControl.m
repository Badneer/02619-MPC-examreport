
%% ---------------------------------------------------------
% Problem 3 - P, PI and PID Control
% Modified Four-Tank System
%
% Input-output pairing:
% F1 -> h1
% F2 -> h2
%
% Experiments:
% 1. Initial-condition recovery
% 2. Disturbance rejection
%
% Constant references in both experiments
% ---------------------------------------------------------

clear; clc; close all;

%% ---------------------------------------------------------
% Add Common folder to MATLAB path
% ---------------------------------------------------------

scriptFolder = fileparts(mfilename('fullpath'));
parentFolder = fileparts(scriptFolder);

addpath(fullfile(parentFolder,'Common'));


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
% Nominal operating conditions
% ---------------------------------------------------------

us = [300;300];        % Nominal pump flows [cm^3/s]
ds = [20;20];          % Nominal disturbances [cm^3/s]

F1s = us(1);
F2s = us(2);

% Steady-state outlet flows
q3s = (1-gamma(2))*F2s + ds(1);
q4s = (1-gamma(1))*F1s + ds(2);

q1s = gamma(1)*F1s + q3s;
q2s = gamma(2)*F2s + q4s;

qs = [q1s;q2s;q3s;q4s];

% Steady-state levels
hs = (qs./a).^2/(2*g);

% Steady-state masses
xs = rho*A.*hs;

% Nominal controlled outputs
zs = FourTankOutput(xs,p);

%% ---------------------------------------------------------
% Steady-state gain matrix and RGA
% ---------------------------------------------------------

G0 = [
    q1s*gamma(1)/(a(1)^2*g), ...
    q1s*(1-gamma(2))/(a(1)^2*g);

    q2s*(1-gamma(1))/(a(2)^2*g), ...
    q2s*gamma(2)/(a(2)^2*g)
];

RGA = G0.*inv(G0)';

disp('Steady-state gain matrix:')
disp(G0)

disp('Relative Gain Array:')
disp(RGA)

%% ---------------------------------------------------------
% Simulation settings
% ---------------------------------------------------------

Ts = 4;                % Controller sampling time [s]
tf = 60*60;            % Simulation duration [s]

T = (0:Ts:tf)';
N = length(T);

% Pump limits [cm^3/s]
umin = [0;0];
umax = [500;500];

%% ---------------------------------------------------------
% Constant reference levels
% ---------------------------------------------------------

R = repmat(zs',N,1);

%% ---------------------------------------------------------
% Controller tuning
% ---------------------------------------------------------

% Initial proportional gains
Kp = 0.5./diag(G0);

% Integral times [s]
Ti = [300;300];

% Derivative times [s]
Td = [15;15];

% Derivative filter time constants
Tf = Td/10;

% Anti-windup tracking time constants
Tt = Ti/2;

% Controller types
controllerNames = {'P','PI','PID'};
nControllers = numel(controllerNames);

%% ---------------------------------------------------------
% Measurement noise
% ---------------------------------------------------------

% Set to zero for the initial stability comparison
sigma_y = [0;0];

% For noisy measurements, use:
% sigma_y = [0.1;0.1];

Rvv = diag(sigma_y.^2);

%% ---------------------------------------------------------
% Experiment definitions
% ---------------------------------------------------------

experimentNames = { ...
    'Initial-condition recovery', ...
    'Disturbance rejection'};

nExperiments = numel(experimentNames);

%% =========================================================
% Run both closed-loop experiments
% ==========================================================

for expID = 1:nExperiments

    %% -----------------------------------------------------
    % Define initial conditions and disturbance sequence
    % -----------------------------------------------------

    D = repmat(ds',N,1);

    if expID == 1

        % Experiment 1:
        % Constant disturbances, displaced initial levels

        h0 = hs;

        % Increase initial levels of Tanks 1 and 2
        h0(1) = hs(1)+2;
        h0(2) = hs(2)+2;

        x0 = rho*A.*h0;

    else

        % Experiment 2:
        % Nominal initial state, step disturbance in F3

        x0 = xs;

        % F3 increases by 10 cm^3/s
        % between 5 and 15 minutes

        for k = 1:N

            if T(k) >= 5*60 && T(k) < 15*60
                D(k,1) = ds(1)+10;
            end

        end

    end

    %% -----------------------------------------------------
    % Storage
    % -----------------------------------------------------

    Xall = zeros(N,4,nControllers);
    Yall = zeros(N,2,nControllers);
    Zall = zeros(N,2,nControllers);
    Uall = zeros(N,2,nControllers);

    %% -----------------------------------------------------
    % Reproducible measurement noise
    % Use the same noise for all controller comparisons
    % -----------------------------------------------------

    rng(100+expID);

    V = randn(N,2)*diag(sigma_y);

    %% =====================================================
    % Closed-loop simulation: P, PI and PID
    % =====================================================

    for c = 1:nControllers

        controllerType = controllerNames{c};

        % Initial physical state
        x = x0;

        % Integral controller state
        I = zeros(2,1);

        % Filtered measured-level derivative
        dyFiltered = zeros(2,1);

        % Previous sensor measurements
        yPrevious = FourTankSensor(x,p) + V(1,:)';

        for k = 1:N

            tk = T(k);

            %% ---------------------------------------------
            % Measurements and controlled outputs
            % ----------------------------------------------

            zk = FourTankOutput(x,p);

            yk = FourTankSensor(x,p) + V(k,:)';

            % Constant reference
            rk = R(k,:)';

            % Control error
            ek = rk-yk;

            %% ---------------------------------------------
            % Filtered derivative of measurement
            % ----------------------------------------------
if k > 1
    dyFiltered = (yk-yPrevious)/Ts;
else
    dyFiltered = zeros(2,1);
end

            %% ---------------------------------------------
            % Controller action
            % ----------------------------------------------

            % Proportional action
            Pterm = Kp.*ek;

            % Integral and derivative action
            Iterm = zeros(2,1);
            Dterm = zeros(2,1);

            switch controllerType

                case 'P'

                    % P controller

                case 'PI'

                    Iterm = I;

                case 'PID'

                    Iterm = I;

                    % Derivative on measurement
                    Dterm = -Kp.*Td.*dyFiltered;

            end

            % Unsaturated command
            uRaw = us + Pterm + Iterm + Dterm;

            % Apply pump limits
            uk = min(max(uRaw,umin),umax);

            %% ---------------------------------------------
            % Store simulation results
            % ----------------------------------------------

            Xall(k,:,c) = x';
            Yall(k,:,c) = yk';
            Zall(k,:,c) = zk';
            Uall(k,:,c) = uk';

            %% ---------------------------------------------
            % Integrate nonlinear plant
            % ----------------------------------------------

            if k < N

                % Disturbance during current interval
                dk = D(k,:)';

                [~,Xseg] = ode15s( ...
                    @(tt,xx) FourTankSystem( ...
                        tt,xx,uk,dk,p), ...
                    [tk T(k+1)],x);

                % Update physical state
                x = Xseg(end,:)';

                %% -----------------------------------------
                % Integral action and anti-windup
                % ------------------------------------------

                if strcmp(controllerType,'PI') || ...
                   strcmp(controllerType,'PID')

                    Ki = Kp./Ti;

                    % Back-calculation anti-windup
                    I = I + Ts*( ...
                        Ki.*ek + (uk-uRaw)./Tt);

                end

            end

            % Update previous measurement
            yPrevious = yk;

        end

    end

    %% =====================================================
    % Figure 1 - Closed-loop tank levels
    % =====================================================

    figure('Name', ...
        [experimentNames{expID} ' - Tank Levels']);

    for j = 1:2

        subplot(2,1,j)
        hold on

        stairs(T/60,R(:,j),'k--', ...
            'LineWidth',1.6, ...
            'DisplayName','Constant reference');

        for c = 1:nControllers

            plot(T/60,Zall(:,j,c), ...
                'LineWidth',1.5, ...
                'DisplayName',controllerNames{c});

        end

        xlabel('Time [min]')
        ylabel(sprintf('h_%d [cm]',j))

        title(sprintf('Tank %d - Controlled Level',j))

        legend('Location','best')
        grid on
        hold off

    end

    sgtitle([experimentNames{expID} ...
        ' - Closed-Loop Responses'])

    %% =====================================================
    % Figure 2 - Manipulated pump flows
    % =====================================================

    figure('Name', ...
        [experimentNames{expID} ' - Pump Flows']);

    for j = 1:2

        subplot(2,1,j)
        hold on

        for c = 1:nControllers

            stairs(T/60,Uall(:,j,c), ...
                'LineWidth',1.5, ...
                'DisplayName',controllerNames{c});

        end

        yline(us(j),'k--', ...
            'LineWidth',1.2, ...
            'DisplayName','Nominal pump flow');

        xlabel('Time [min]')
        ylabel(sprintf('F_%d [cm^3/s]',j))

        title(sprintf('Pump %d',j))

        legend('Location','best')
        grid on
        hold off

    end

    sgtitle([experimentNames{expID} ...
        ' - Manipulated Inputs'])

    %% =====================================================
    % Figure 3 - Tracking errors
    % =====================================================

    figure('Name', ...
        [experimentNames{expID} ' - Tracking Errors']);

    for j = 1:2

        subplot(2,1,j)
        hold on

        for c = 1:nControllers

            errorSignal = R(:,j)-Zall(:,j,c);

            plot(T/60,errorSignal, ...
                'LineWidth',1.5, ...
                'DisplayName',controllerNames{c});

        end

        yline(0,'k--', ...
            'DisplayName','Zero error');

        xlabel('Time [min]')
        ylabel(sprintf('e_%d [cm]',j))

        title(sprintf('Tank %d - Tracking Error',j))

        legend('Location','best')
        grid on
        hold off

    end

    sgtitle([experimentNames{expID} ...
        ' - Tracking Errors'])

    %% =====================================================
    % Figure 4 - Applied disturbance flows
    % =====================================================

    figure('Name', ...
        [experimentNames{expID} ' - Disturbances']);

    for j = 1:2

        subplot(2,1,j)

        stairs(T/60,D(:,j), ...
            'LineWidth',1.5)

        hold on

        yline(ds(j),'k--', ...
            'LineWidth',1.2)

        xlabel('Time [min]')
        ylabel(sprintf('F_%d [cm^3/s]',j+2))

        title(sprintf('Disturbance F%d',j+2))

        legend('Applied disturbance','Nominal', ...
            'Location','best')

        grid on

    end

    sgtitle([experimentNames{expID} ...
        ' - Disturbance Inputs'])

    %% =====================================================
    % Performance metrics
    % =====================================================

    fprintf('\n========================================\n');
    fprintf('Experiment %d: %s\n', ...
        expID,experimentNames{expID});
    fprintf('========================================\n');

    for c = 1:nControllers

        E = R-Zall(:,:,c);

        % Integral of Absolute Error
        IAE = trapz(T,abs(E));

        % Root Mean Square Error
        RMSE = sqrt(mean(E.^2,1));

        % Final tracking error
        finalError = E(end,:);

        fprintf('\n%s Controller\n',controllerNames{c});

        fprintf('Tank 1 IAE = %.3f cm*s\n',IAE(1));
        fprintf('Tank 2 IAE = %.3f cm*s\n',IAE(2));

        fprintf('Tank 1 RMSE = %.4f cm\n',RMSE(1));
        fprintf('Tank 2 RMSE = %.4f cm\n',RMSE(2));

        fprintf('Final error Tank 1 = %.4f cm\n', ...
            finalError(1));

        fprintf('Final error Tank 2 = %.4f cm\n', ...
            finalError(2));

    end

end

%% ---------------------------------------------------------
% Display controller parameters
% ---------------------------------------------------------

fprintf('\n----------------------------------------\n');
fprintf('Controller tuning parameters\n');
fprintf('----------------------------------------\n');

disp('Proportional gains Kp:')
disp(Kp)

disp('Integral times Ti [s]:')
disp(Ti)

disp('Derivative times Td [s]:')
disp(Td)

disp('Measurement covariance Rvv:')
disp(Rvv)
