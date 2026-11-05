
%% Problem 4.1 and 4.3
% Deterministic nonlinear step responses
% Modified Four-Tank System

clear; clc; close all;

%% Add Common folder from Problem 2

scriptFolder = fileparts(mfilename('fullpath'));

% Move to the project root
projectFolder = fileparts(scriptFolder);

% Locate Common inside Problem 2
commonFolder = fullfile(projectFolder, 'Problem 2', 'Common');

if ~isfolder(commonFolder)
    error('Common folder not found: %s', commonFolder);
end

addpath(commonFolder);

%% --------------------------------------------------------
%  Physical parameters
% ---------------------------------------------------------

a = 1.2272*ones(4,1);       % Outlet areas [cm^2]
A = 380.1327*ones(4,1);     % Tank areas [cm^2]

gamma = [0.6;0.7];

g = 981;                     % Gravity [cm/s^2]
rho = 1.00;                  % Density [g/cm^3]

p = [a;A;gamma;g;rho];

%% --------------------------------------------------------
%  Nominal operating point
% ---------------------------------------------------------

us = [300;300];              % Pump flows [cm^3/s]
ds = [20;20];                % Disturbance flows [cm^3/s]

% Steady-state outlet flows

q3s = (1-gamma(2))*us(2) + ds(1);
q4s = (1-gamma(1))*us(1) + ds(2);

q1s = gamma(1)*us(1) + q3s;
q2s = gamma(2)*us(2) + q4s;

qs = [q1s;q2s;q3s;q4s];

% Steady-state tank levels
hs = (qs./a).^2/(2*g);

% Steady-state masses
xs = rho*A.*hs;

% Nominal controlled outputs
ys = FourTankOutput(xs,p);

fprintf('\nNominal steady-state levels:\n');
fprintf('h1 = %.4f cm\n',ys(1));
fprintf('h2 = %.4f cm\n',ys(2));

%% --------------------------------------------------------
%  Simulation settings
% ---------------------------------------------------------

Ts = 4;                       % Sampling time [s]
tf = 60*60;                   % Final time: 60 min [s]
tStep = 5*60;                 % Step at 5 min [s]

T = (0:Ts:tf)';
Nt = length(T);

stepSizes = [0.10 0.25 0.50];

nSteps = length(stepSizes);
nInputs = 2;
nOutputs = 2;

% Storage:
% Y(time,output,step size,input)
% U(time,input,step size,stepped input)

Y = zeros(Nt,nOutputs,nSteps,nInputs);
U = zeros(Nt,nInputs,nSteps,nInputs);

%% --------------------------------------------------------
%  Problem 4.1: Nonlinear step experiments
% ---------------------------------------------------------

for j = 1:nInputs

    for k = 1:nSteps

        deltaU = stepSizes(k)*us(j);

        % Start from nominal steady state
        x = xs;

        Y(1,:,k,j) = FourTankOutput(x,p)';
        U(1,:,k,j) = us';

        for n = 1:Nt-1

            % Piecewise-constant manipulated input
            uk = us;

            if T(n) >= tStep
                uk(j) = us(j) + deltaU;
            end

            % Nonlinear deterministic plant
            [~,Xseg] = ode15s( ...
                @(t,x) FourTankSystem(t,x,uk,ds,p), ...
                [T(n) T(n+1)],x);

            x = Xseg(end,:)';

            % Controlled output
            Y(n+1,:,k,j) = FourTankOutput(x,p)';

            % Manipulated input history
            U(n+1,:,k,j) = uk';

        end

    end

end

%% --------------------------------------------------------
%  Plot 1: Responses to F1
%  Plot 2: Responses to F2
% ---------------------------------------------------------

for j = 1:nInputs

    figure('Name',sprintf('Nonlinear step responses - F%d',j), ...
        'Color','w');

    for i = 1:nOutputs

        subplot(2,1,i);
        hold on;
        grid on;
        box on;

        for k = 1:nSteps

            plot(T/60,Y(:,i,k,j), ...
                'LineWidth',1.6);

        end

        yline(ys(i),'k--','Nominal');

        xlabel('Time [min]');
        ylabel(sprintf('h_%d [cm]',i));

        title(sprintf('Response of Tank %d to F%d',i,j));

        legend('+10%','+25%','+50%', ...
            'Nominal','Location','best');

    end

end

%% --------------------------------------------------------
%  Problem 4.3: Normalized step responses
%
%  S_ij(t) = (y_i(t) - y_i,s) / Delta u_j
%
%  Units: cm / (cm^3/s)
% ---------------------------------------------------------

Ynorm = zeros(size(Y));

for j = 1:nInputs

    for k = 1:nSteps

        deltaU = stepSizes(k)*us(j);

        for i = 1:nOutputs

            Ynorm(:,i,k,j) = ...
                (Y(:,i,k,j)-ys(i))/deltaU;

        end

    end

end

%% --------------------------------------------------------
%  Plot 3: Normalized responses to F1
%  Plot 4: Normalized responses to F2
% ---------------------------------------------------------

for j = 1:nInputs

    figure('Name',sprintf('Normalized step responses - F%d',j), ...
        'Color','w');

    for i = 1:nOutputs

        subplot(2,1,i);
        hold on;
        grid on;
        box on;

        for k = 1:nSteps

            plot(T/60,Ynorm(:,i,k,j), ...
                'LineWidth',1.6);

        end

        xlabel('Time [min]');

        ylabel(sprintf('Delta h_%d / Delta F_%d',i,j));

        title(sprintf( ...
            'Normalized response: F%d to h%d',j,i));

        legend('+10%','+25%','+50%', ...
            'Location','best');

    end

end

%% --------------------------------------------------------
%  Steady-state gains from normalized steps
% ---------------------------------------------------------

fprintf('\n-----------------------------------------\n');
fprintf('STEADY-STATE GAINS FROM STEP EXPERIMENTS\n');
fprintf('-----------------------------------------\n');

for k = 1:nSteps

    Gss = zeros(nOutputs,nInputs);

    for j = 1:nInputs
        for i = 1:nOutputs

            Gss(i,j) = Ynorm(end,i,k,j);

        end
    end

    fprintf('\nStep size: +%.0f%%\n',100*stepSizes(k));
    disp(Gss);

end

fprintf('\nSimulation completed.\n');
%% Save deterministic normalized responses

scriptFolder = fileparts(mfilename('fullpath'));

save(fullfile(scriptFolder,'StepData41.mat'), ...
    'T','tStep','stepSizes','Ynorm');