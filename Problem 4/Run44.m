
%% Problem 4.4
% Transfer Function Identification
% Modified Four-Tank System
% First-order direct paths
% Second-order cross-coupling paths

clear; clc; close all;

%% --------------------------------------------------------
% Load results from Problems 4.1 and 4.2
% ---------------------------------------------------------

scriptFolder = fileparts(mfilename('fullpath'));

det = load(fullfile(scriptFolder,'StepData41.mat'));
noi = load(fullfile(scriptFolder,'StepData42.mat'));

T = det.T(:);
tStep = det.tStep;
stepSizes = det.stepSizes;

Ydet = det.Ynorm;
Ynoisy = noi.YnormNoisy;

noiseNames = cellstr(noi.noiseNames);

nSteps = length(stepSizes);
nNoise = length(noiseNames);

%% --------------------------------------------------------
% Identification settings
% ---------------------------------------------------------

nInputs = 2;
nOutputs = 2;

steadyWindow = 5*60;

idxSS = T >= (T(end)-steadyWindow);
idxPre = T < tStep;
idxPost = T >= tStep;

tRel = T(idxPost)-tStep;

% Dimensions: output, input, step, scenario
K = zeros(2,2,nSteps,nNoise+1);
tau = zeros(2,2,nSteps,nNoise+1);

%% --------------------------------------------------------
% Identify gains and 63.2% characteristic times
% ---------------------------------------------------------

for scenario = 1:nNoise+1

    for k = 1:nSteps

        for j = 1:nInputs

            for i = 1:nOutputs

                if scenario == 1
                    y = Ydet(:,i,k,j);
                else
                    y = Ynoisy(:,i,k,j,scenario-1);
                end

                y = y(:);

                % Initial response level
                y0 = mean(y(idxPre));

                % Steady-state change
                ySS = mean(y(idxSS));

                Kij = ySS-y0;

                K(i,j,k,scenario) = Kij;

                % 63.2% response level
                target = y0+0.632*Kij;

                yPost = y(idxPost);
                tPost = T(idxPost);

                if Kij >= 0
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
% Display identified parameters
% ---------------------------------------------------------

scenarioNames = [{'Deterministic'},noiseNames(:)'];

for scenario = 1:nNoise+1

    fprintf('\n========================================\n');
    fprintf('%s\n',scenarioNames{scenario});
    fprintf('========================================\n');

    for k = 1:nSteps

        fprintf('\nStep size: +%.0f%%\n', ...
                100*stepSizes(k));

        fprintf('\nSteady-state gain matrix K:\n');
        disp(K(:,:,k,scenario));

        fprintf('63.2%% characteristic times [s]:\n');
        disp(tau(:,:,k,scenario));

    end
end

%% --------------------------------------------------------
% Reference model: Deterministic 10% step
% ---------------------------------------------------------

kRef = 1;
scenarioRef = 1;

Kref = K(:,:,kRef,scenarioRef);
tauRef = tau(:,:,kRef,scenarioRef);

% First-order:
% G11 and G22
%
% Second-order:
% G12 and G21

tauSecond = NaN(2,2,2);

%% --------------------------------------------------------
% Second-order parameter fitting
% ---------------------------------------------------------

% The second-order step response is:
%
% y(t) = K * [1 -
% (tauA*exp(-t/tauA)-tauB*exp(-t/tauB))
% /(tauA-tauB)]
%
% Fit tauA and tauB by minimizing squared errors.
% K is fixed to the identified steady-state gain.

crossPaths = [1 2;
              2 1];

options = optimset('Display','off', ...
                   'TolX',1e-8, ...
                   'TolFun',1e-10, ...
                   'MaxIter',3000, ...
                   'MaxFunEvals',6000);

for path = 1:2

    i = crossPaths(path,1);
    j = crossPaths(path,2);

    y = Ydet(:,i,kRef,j);
    y = y(:);

    y0 = mean(y(idxPre));

    yFit = y(idxPost)-y0;

    Kij = Kref(i,j);

    % Initial guesses based on the 63.2% time
    tau63 = tauRef(i,j);

    if ~isfinite(tau63) || tau63 <= 0
        tau63 = 100;
    end

    tauGuess = tau63/2.15;

    p0 = log([0.7*tauGuess;
              1.3*tauGuess]);

    % Optimization in logarithmic coordinates
    % ensures positive time constants.

    objective = @(p) sum( ...
        (secondOrderStep(exp(p),tRel,Kij) ...
        -yFit).^2) ...
        + 1e4*sum(max(0,abs(p)-log(1e5)).^2);

    pOpt = fminsearch(objective,p0,options);

    tauFit = sort(exp(pOpt));

    tauSecond(i,j,1) = tauFit(1);
    tauSecond(i,j,2) = tauFit(2);

    fprintf('\n========================================\n');
    fprintf('Second-order fitting: G%d%d\n',i,j);
    fprintf('========================================\n');

    fprintf('K = %.6f\n',Kij);
    fprintf('tauA = %.4f s\n',tauFit(1));
    fprintf('tauB = %.4f s\n',tauFit(2));

end

%% --------------------------------------------------------
% Construct transfer function matrix
% ---------------------------------------------------------

s = tf('s');

% First-order direct transfer functions

G11 = Kref(1,1)/(tauRef(1,1)*s+1);

G22 = Kref(2,2)/(tauRef(2,2)*s+1);

% Second-order cross transfer functions

tau12a = tauSecond(1,2,1);
tau12b = tauSecond(1,2,2);

tau21a = tauSecond(2,1,1);
tau21b = tauSecond(2,1,2);

G12 = Kref(1,2)/ ...
    ((tau12a*s+1)*(tau12b*s+1));

G21 = Kref(2,1)/ ...
    ((tau21a*s+1)*(tau21b*s+1));

G = [G11 G12;
     G21 G22];

fprintf('\n========================================\n');
fprintf('IDENTIFIED TRANSFER FUNCTION MATRIX\n');
fprintf('Deterministic +10%% step\n');
fprintf('========================================\n');

disp(G);

%% --------------------------------------------------------
% Model validation against nonlinear step responses
% ---------------------------------------------------------

figure('Name','Identified Model Validation', ...
       'Color','w');

RMSE = zeros(2,2);
FitPercent = zeros(2,2);

for i = 1:2

    for j = 1:2

        subplot(2,2,(i-1)*2+j);
        hold on; grid on; box on;

        % Nonlinear normalized response
        y = Ydet(:,i,kRef,j);
        y = y(:);

        y0 = mean(y(idxPre));
        yActual = y(idxPost)-y0;

        % Identified model step response
        Gij = G(i,j);

        yModel = step(Gij,tRel);
        yModel = yModel(:);

        % Error metrics
        errorVector = yActual-yModel;

        RMSE(i,j) = sqrt(mean(errorVector.^2));

        denominator = norm(yActual-mean(yActual));

        if denominator > 0
            FitPercent(i,j) = 100* ...
                (1-norm(errorVector)/denominator);
        else
            FitPercent(i,j) = NaN;
        end

        plot(tRel/60,yActual, ...
            'b','LineWidth',1.5);

        plot(tRel/60,yModel, ...
            'r--','LineWidth',1.5);

        xlabel('Time after step [min]');
        ylabel('Normalized response');

        if i == j
            title(sprintf('G%d%d - First Order',i,j));
        else
            title(sprintf('G%d%d - Second Order',i,j));
        end

        legend('Nonlinear model', ...
               'Identified model', ...
               'Location','best');

    end
end

sgtitle('Nonlinear vs Identified Linear Model');

%% --------------------------------------------------------
% Display fitting accuracy
% ---------------------------------------------------------

fprintf('\n========================================\n');
fprintf('MODEL FITTING ACCURACY\n');
fprintf('========================================\n');

fprintf('\nRMSE matrix:\n');
disp(RMSE);

fprintf('\nModel fit percentage:\n');
disp(FitPercent);

%% --------------------------------------------------------
% Save identified model
% ---------------------------------------------------------

save(fullfile(scriptFolder,'IdentifiedModel44.mat'), ...
    'K','tau','Kref','tauRef', ...
    'tauSecond','G','RMSE','FitPercent');

fprintf('\nProblem 4.4 completed.\n');

%% --------------------------------------------------------
% Local function: Second-order step response
% ---------------------------------------------------------

function y = secondOrderStep(tauPair,t,K)

    tauA = tauPair(1);
    tauB = tauPair(2);

    t = t(:);

    % Handle nearly identical time constants
    if abs(tauA-tauB) < ...
            1e-5*max(tauA,tauB)

        tauMean = (tauA+tauB)/2;

        y = K*(1-(1+t/tauMean) ...
            .*exp(-t/tauMean));

    else

        y = K*(1- ...
            (tauA*exp(-t/tauA) ...
            -tauB*exp(-t/tauB)) ...
            /(tauA-tauB));

    end

end
