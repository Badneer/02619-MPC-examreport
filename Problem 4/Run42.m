%% Problem 4.2 and 4.3
% Deterministic nonlinear step responses
% Modified Four-Tank System
% Three levels of noise
% 1) low noise -> sigma_F3 = sigma_F4 = 1 cm3/s (Process Noise)-
% sigma_y = 0.02 cm (Measurement Noise)
% 2) medium noise -> sigma_F3 = sigma_F4 = 5 cm3/s (Process Noise)-
% sigma_y = 0.1 cm (Measurement Noise)
% 3) high noise -> sigma_F3 = sigma_F4 = 10 cm3/s (Process Noise)-
% sigma_y = 0.2 cm (Measurement Noise)

clear; clc; close all;

%% --------------------------------------------------------
% Physical parameters
% ---------------------------------------------------------

a = 1.2272*ones(4,1);
A = 380.1327*ones(4,1);

gamma = [0.6;0.7];

g = 981;
rho = 1.00;

p = [a;A;gamma;g;rho];

%% --------------------------------------------------------
% Nominal operating point
% ---------------------------------------------------------

us = [300;300];        % Pump flows [cm^3/s]
ds = [20;20];          % Mean disturbances [cm^3/s]

q3s = (1-gamma(2))*us(2) + ds(1);
q4s = (1-gamma(1))*us(1) + ds(2);

q1s = gamma(1)*us(1) + q3s;
q2s = gamma(2)*us(2) + q4s;

qs = [q1s;q2s;q3s;q4s];

hs = (qs./a).^2/(2*g);
xs = rho*A.*hs;

ys = FourTankOutput(xs,p);

fprintf('\nNominal steady-state levels:\n');
fprintf('h1 = %.4f cm\n',ys(1));
fprintf('h2 = %.4f cm\n',ys(2));

%% --------------------------------------------------------
% Simulation settings
% ---------------------------------------------------------

Ts = 4;                % Sampling time [s]
tf = 60*60;            % Final time [s]
tStep = 5*60;          % Step at 5 min [s]

T = (0:Ts:tf)';
Nt = length(T);

stepSizes = [0.10 0.25 0.50];

nSteps = length(stepSizes);
nInputs = 2;
nOutputs = 2;

%% --------------------------------------------------------
% Noise scenarios
% ---------------------------------------------------------

noiseNames = {'Low','Medium','High'};

% Standard deviations
sigmaF = [1 5 10];          % F3 and F4 [cm^3/s]
sigmaY = [0.02 0.1 0.2];    % Measurements [cm]

nNoise = length(noiseNames);

rng(100);

% Generate independent noise for each scenario
% Process noise updated every Ts seconds

processNoise = zeros(Nt-1,2,nNoise);
measurementNoise = zeros(Nt,2,nNoise);

for s = 1:nNoise

    processNoise(:,:,s) = sigmaF(s)*randn(Nt-1,2);

    measurementNoise(:,:,s) = sigmaY(s)*randn(Nt,2);

end

%% --------------------------------------------------------
% Storage
% Ydet(time,output,step,input)
% Ynoisy(time,output,step,input,noise scenario)
% ---------------------------------------------------------

Ydet = zeros(Nt,nOutputs,nSteps,nInputs);

Ynoisy = zeros(Nt,nOutputs,nSteps,nInputs,nNoise);

U = zeros(Nt,nInputs,nSteps,nInputs);

%% ========================================================
% Problem 4.1 - Nonlinear step experiments
% =========================================================

for j = 1:nInputs

    for k = 1:nSteps

        deltaU = stepSizes(k)*us(j);

        %% Deterministic simulation

        x = xs;
        Ydet(1,:,k,j) = FourTankOutput(x,p)';

        for n = 1:Nt-1

            uk = us;

            if T(n) >= tStep
                uk(j) = us(j) + deltaU;
            end

            [~,Xseg] = ode15s( ...
                @(t,x) FourTankSystem(t,x,uk,ds,p), ...
                [T(n) T(n+1)],x);

            x = Xseg(end,:)';

            Ydet(n+1,:,k,j) = FourTankOutput(x,p)';

        end

        % Manipulated input history
        U(:,:,k,j) = repmat(us',Nt,1);
        U(T>=tStep,j,k,j) = us(j) + deltaU;

        %% Stochastic simulations

        for s = 1:nNoise

            x = xs;

            % Initial noisy measurement
            Ynoisy(1,:,k,j,s) = ...
                FourTankOutput(x,p)' + measurementNoise(1,:,s);

            for n = 1:Nt-1

                % Manipulated input
                uk = us;

                if T(n) >= tStep
                    uk(j) = us(j) + deltaU;
                end

                % Process disturbances F3 and F4
                dk = ds + processNoise(n,:,s)';

                % Prevent negative inlet flows
                dk = max(dk,0);

                % Nonlinear plant
                [~,Xseg] = ode15s( ...
                    @(t,x) FourTankSystem(t,x,uk,dk,p), ...
                    [T(n) T(n+1)],x);

                x = Xseg(end,:)';

                % Measurement noise
                Ynoisy(n+1,:,k,j,s) = ...
                    FourTankOutput(x,p)' + measurementNoise(n+1,:,s);

            end

        end

    end

end

%% ========================================================
% Problem 4.1 - Plot deterministic step responses
% =========================================================

for j = 1:nInputs

    figure('Name',sprintf('Deterministic responses - F%d',j), ...
        'Color','w');

    for i = 1:nOutputs

        subplot(2,1,i);
        hold on; grid on; box on;

        for k = 1:nSteps

            plot(T/60,Ydet(:,i,k,j),'LineWidth',1.6);

        end

        yline(ys(i),'k--','Nominal');

        xlabel('Time [min]');
        ylabel(sprintf('h_%d [cm]',i));

        title(sprintf('Deterministic response: F%d to h%d',j,i));

        legend('+10%','+25%','+50%','Nominal', ...
            'Location','best');

    end

end

%% ========================================================
% Problem 4.1 - Plot noisy step responses
% One figure per input and noise scenario
% =========================================================

for s = 1:nNoise

    for j = 1:nInputs

        figure('Name',sprintf('%s noise - F%d',noiseNames{s},j), ...
            'Color','w');

        for i = 1:nOutputs

            subplot(2,1,i);
            hold on; grid on; box on;

            for k = 1:nSteps

                plot(T/60,Ynoisy(:,i,k,j,s), ...
                    'LineWidth',1.0);

            end

            yline(ys(i),'k--','Nominal');

            xlabel('Time [min]');
            ylabel(sprintf('h_%d [cm]',i));

            title(sprintf('%s noise: F%d to h%d', ...
                noiseNames{s},j,i));

            legend('+10%','+25%','+50%','Nominal', ...
                'Location','best');

        end

    end

end

%% ========================================================
% Problem 4.3 - Normalized step responses
%
% S_ij(t) = (y_i(t)-y_i,s)/Delta u_j
% =========================================================

YnormDet = zeros(size(Ydet));
YnormNoisy = zeros(size(Ynoisy));

for j = 1:nInputs

    for k = 1:nSteps

        deltaU = stepSizes(k)*us(j);

        for i = 1:nOutputs

            % Deterministic
            YnormDet(:,i,k,j) = ...
                (Ydet(:,i,k,j)-ys(i))/deltaU;

            % Noisy
            for s = 1:nNoise

                YnormNoisy(:,i,k,j,s) = ...
                    (Ynoisy(:,i,k,j,s)-ys(i))/deltaU;

            end

        end

    end

end

%% ========================================================
% Problem 4.3 - Plot deterministic normalized responses
% =========================================================

for j = 1:nInputs

    figure('Name',sprintf('Normalized deterministic - F%d',j), ...
        'Color','w');

    for i = 1:nOutputs

        subplot(2,1,i);
        hold on; grid on; box on;

        for k = 1:nSteps

            plot(T/60,YnormDet(:,i,k,j), ...
                'LineWidth',1.6);

        end

        xlabel('Time [min]');
        ylabel(sprintf('Delta h_%d / Delta F_%d',i,j));

        title(sprintf('Normalized deterministic: F%d to h%d',j,i));

        legend('+10%','+25%','+50%', ...
            'Location','best');

    end

end

%% ========================================================
% Problem 4.3 - Plot noisy normalized responses
% =========================================================

for s = 1:nNoise

    for j = 1:nInputs

        figure('Name',sprintf('Normalized %s noise - F%d', ...
            noiseNames{s},j),'Color','w');

        for i = 1:nOutputs

            subplot(2,1,i);
            hold on; grid on; box on;

            for k = 1:nSteps

                plot(T/60,YnormNoisy(:,i,k,j,s), ...
                    'LineWidth',1.0);

            end

            xlabel('Time [min]');
            ylabel(sprintf('Delta h_%d / Delta F_%d',i,j));

            title(sprintf('Normalized %s noise: F%d to h%d', ...
                noiseNames{s},j,i));

            legend('+10%','+25%','+50%', ...
                'Location','best');

        end

    end

end

%% ========================================================
% Steady-state gains - Deterministic reference
% =========================================================

fprintf('\n-----------------------------------------\n');
fprintf('DETERMINISTIC STEADY-STATE GAINS\n');
fprintf('-----------------------------------------\n');

for k = 1:nSteps

    Gss = zeros(nOutputs,nInputs);

    for j = 1:nInputs
        for i = 1:nOutputs

            Gss(i,j) = YnormDet(end,i,k,j);

        end
    end

    fprintf('\nStep size: +%.0f%%\n',100*stepSizes(k));
    disp(Gss);

end

%% ========================================================
% Noise scenario information
% =========================================================

fprintf('\n-----------------------------------------\n');
fprintf('NOISE SCENARIOS\n');
fprintf('-----------------------------------------\n');

for s = 1:nNoise

    fprintf('%s: sigmaF = %.2f cm^3/s, sigmaY = %.2f cm\n', ...
        noiseNames{s},sigmaF(s),sigmaY(s));

end

fprintf('\nSimulation completed.\n');

%% Save noisy normalized responses

scriptFolder = fileparts(mfilename('fullpath'));

save(fullfile(scriptFolder,'StepData42.mat'), ...
    'T','tStep','stepSizes','YnormNoisy','noiseNames');