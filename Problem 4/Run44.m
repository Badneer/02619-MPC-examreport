
%% Problem 4.4
% Transfer Function Identification
% Modified Four-Tank System
% Deterministic and Stochastic Step Responses

clear; clc; close all;

%% --------------------------------------------------------
% Load results from Problems 4.1 and 4.2
% ---------------------------------------------------------

scriptFolder = fileparts(mfilename('fullpath'));

det = load(fullfile(scriptFolder,'StepData41.mat'));
noi = load(fullfile(scriptFolder,'StepData42.mat'));

T = det.T;
tStep = det.tStep;
stepSizes = det.stepSizes;

Ydet = det.Ynorm;
Ynoisy = noi.YnormNoisy;

noiseNames = noi.noiseNames;

nSteps = length(stepSizes);
nNoise = length(noiseNames);

%% --------------------------------------------------------
% Identification settings
% ---------------------------------------------------------

nInputs = 2;
nOutputs = 2;

% Average over the last 5 minutes for steady-state estimate
steadyWindow = 5*60;

idxSS = T >= (T(end)-steadyWindow);

% Store identified parameters
% Dimensions: output, input, step size, scenario
% Scenario 1 = deterministic
% Scenario 2 = low noise
% Scenario 3 = medium noise
% Scenario 4 = high noise

K = zeros(2,2,nSteps,nNoise+1);
tau = zeros(2,2,nSteps,nNoise+1);

%% --------------------------------------------------------
% Identify parameters
% ---------------------------------------------------------

for scenario = 1:nNoise+1

    for k = 1:nSteps

        for j = 1:nInputs

            for i = 1:nOutputs

                % Select normalized response
                if scenario == 1
                    y = Ydet(:,i,k,j);
                else
                    y = Ynoisy(:,i,k,j,scenario-1);
                end

                y = y(:);

                % Initial level before step
                idx0 = find(T < tStep);

                y0 = mean(y(idx0));

                % Steady-state gain
                Kij = mean(y(idxSS));

                K(i,j,k,scenario) = Kij;

                % Target at 63.2% of total change
                target = y0 + 0.632*(Kij-y0);

                % Find crossing after step
                idxPost = find(T >= tStep);

                yPost = y(idxPost);
                tPost = T(idxPost);

                if Kij >= y0
                    idx = find(yPost >= target,1,'first');
                else
                    idx = find(yPost <= target,1,'first');
                end

                if isempty(idx)
                    tau(i,j,k,scenario) = NaN;
                    continue;
                end

                % Linear interpolation
                if idx == 1
                    t63 = tPost(1);
                else

                    t1 = tPost(idx-1);
                    t2 = tPost(idx);

                    y1 = yPost(idx-1);
                    y2 = yPost(idx);

                    if abs(y2-y1) < 1e-12
                        t63 = t2;
                    else
                        t63 = t1 + ...
                            (target-y1)*(t2-t1)/(y2-y1);
                    end

                end

                tau(i,j,k,scenario) = t63-tStep;

            end

        end

    end

end

%% --------------------------------------------------------
% Display identified gains and time constants
% ---------------------------------------------------------

scenarioNames = [{'Deterministic'},noiseNames];

for scenario = 1:nNoise+1

    fprintf('\n========================================\n');
    fprintf('%s\n',scenarioNames{scenario});
    fprintf('========================================\n');

    for k = 1:nSteps

        fprintf('\nStep size: +%.0f%%\n',100*stepSizes(k));

        fprintf('\nSteady-state gain matrix K:\n');
        disp(K(:,:,k,scenario));

        fprintf('Time constant matrix tau [s]:\n');
        disp(tau(:,:,k,scenario));

    end

end

%% --------------------------------------------------------
% Transfer functions
% Reference: deterministic 10% step
% ---------------------------------------------------------

kRef = 1;
scenarioRef = 1;

Kref = K(:,:,kRef,scenarioRef);
tauRef = tau(:,:,kRef,scenarioRef);

G11 = tf(Kref(1,1),[tauRef(1,1) 1]);
G12 = tf(Kref(1,2),[tauRef(1,2) 1]);

G21 = tf(Kref(2,1),[tauRef(2,1) 1]);
G22 = tf(Kref(2,2),[tauRef(2,2) 1]);

G = [G11 G12;
     G21 G22];

fprintf('\n========================================\n');
fprintf('REFERENCE TRANSFER FUNCTION MATRIX\n');
fprintf('Deterministic +10%% step\n');
fprintf('========================================\n');

disp(G);

%% --------------------------------------------------------
% Compare identified model with deterministic responses
% ---------------------------------------------------------

tRel = T-tStep;
idxPlot = T >= tStep;

figure('Name','Identified Model Validation', ...
       'Color','w');

for i = 1:2

    for j = 1:2

        subplot(2,2,(i-1)*2+j);
        hold on; grid on; box on;

        % Actual normalized response
        yActual = Ydet(idxPlot,i,kRef,j);

        % First-order identified model response
        yModel = Kref(i,j)* ...
            (1-exp(-tRel(idxPlot)/tauRef(i,j)));

        plot(tRel(idxPlot)/60,yActual, ...
            'b','LineWidth',1.5);

        plot(tRel(idxPlot)/60,yModel, ...
            'r--','LineWidth',1.5);

        xlabel('Time after step [min]');
        ylabel('Normalized response');

        title(sprintf('G%d%d',i,j));

        legend('Nonlinear model','Identified model', ...
               'Location','best');

    end

end

%% --------------------------------------------------------
% Save identified model
% ---------------------------------------------------------

save(fullfile(scriptFolder,'IdentifiedModel44.mat'), ...
    'K','tau','Kref','tauRef','G');

fprintf('\nProblem 4.4 completed.\n');
