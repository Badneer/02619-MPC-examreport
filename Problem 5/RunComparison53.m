%% Problem 5.3 - Compare linearization with the identified step models

% Problem 4: Run41.m, Run42.m and Run44.m.
% Uses their saved gains/time constants without performing a new fit.

scriptFolder = fileparts(mfilename('fullpath'));
projectFolder = fileparts(scriptFolder);
problem4Folder = fullfile(projectFolder,'Problem 4');

stepDet = load(fullfile(problem4Folder,'StepData41.mat'));
stepNoisy = load(fullfile(problem4Folder,'StepData42.mat'));
identified = load(fullfile(problem4Folder,'IdentifiedModel44.mat'));

% The identification file uses the order of steps/scenarios in these files.
if ~isequal(stepDet.T,stepNoisy.T) || ...
   ~isequal(stepDet.tStep,stepNoisy.tStep) || ...
   ~isequal(stepDet.stepSizes,stepNoisy.stepSizes)
    error('Problem 4 datasets have inconsistent times or step sizes.');
end
scenarioNames = [{'Deterministic'},stepNoisy.noiseNames];
stepSizes = stepDet.stepSizes;
if size(identified.K,3) ~= numel(stepSizes) || ...
   size(identified.K,4) ~= numel(scenarioNames) || ...
   ~isequal(size(identified.K),size(identified.tau))
    error('Identified parameters do not match the Problem 4 data dimensions.');
end

%% Compute the continuous-time models from Problem 5.1
run(fullfile(scriptFolder,'RunLinearization51.m'));

% Run41.m and Run42.m explicitly use us=[300;300], ds=[20;20]
% and the same physical parameters as Problem 5.1.
% sysPC is the four-state plant at exactly that operating point.

% Its pump transfer functions also equal the pump channels of sysSDE.
% sysDet uses dsDet=[0;0], so it is not a matched operating-point comparison.
pProblem4 = [1.2272*ones(4,1);380.1327*ones(4,1);0.6;0.7;981;1.00];
if ~isequal(us,[300;300]) || ~isequal(dsPC,[20;20]) || ~isequal(p,pProblem4)
    error('Problem 5 parameters/nominal flows do not match Run41.m and Run42.m.');
end

Klinear = dcgain(sysPC(:,1:2));
TauTank = -1./diag(AcPC);

fprintf('\n5.3: matched operating point: F1=F2=300, F3=F4=20 cm^3/s\n');
fprintf('Linearized pump-to-level gains (sysPC and sysSDE):\n');
disp(Klinear)
fprintf('Tank pole time constants [Tau1;Tau2;Tau3;Tau4], seconds:\n');
disp(TauTank)
fprintf('Direct channels: G11 has Tau1; G22 has Tau2.\n');
fprintf('Cross channels: G12 has Tau1 and Tau3; G21 has Tau2 and Tau4.\n');
fprintf('sysDet is excluded from percentage comparisons: its nominal F3=F4=0.\n');

%% Comparable 63.2%% response times for the linearized plant

% Run44 identifies ONE first-order time constant from the first crossing
% of 63.2%% of the total response. The exact cross channels have TWO poles.

% Apply the same crossing rule and sampling grid to the exact linear step
% responses. This does not replace either of their actual pole time constants.

tRelative = stepDet.T(stepDet.T >= stepDet.tStep)-stepDet.tStep;
Slinear = zeros(numel(tRelative),2,2);
t632Linear = NaN(2,2);

for i = 1:2
    for j = 1:2
        tauLower = TauTank(i);
        if i == j
            shape = 1-exp(-tRelative/tauLower);
        else
            tauUpper = TauTank(i+2);
            shape = 1-(tauLower*exp(-tRelative/tauLower) ...
                -tauUpper*exp(-tRelative/tauUpper))/(tauLower-tauUpper);
        end
        Slinear(:,i,j) = Klinear(i,j)*shape;
        target = 0.632*Klinear(i,j);
        idx = find(Slinear(:,i,j) >= target,1,'first');
        if isempty(idx)
            error('Linear response does not reach 63.2%% within the recorded horizon.');
        end
        t1 = tRelative(idx-1);
        t2 = tRelative(idx);
        y1 = Slinear(idx-1,i,j);
        y2 = Slinear(idx,i,j);
        t632Linear(i,j) = t1+(target-y1)*(t2-t1)/(y2-y1);
    end
end

fprintf('\nLinearized 63.2%% response times [s], same rule as Run44:\n');
disp(t632Linear)

%% Compare every recorded step size and noise scenario
% K4: gain identified in Run44; K5: gain of the matched linearization.
% tau4: first-order time constant identified in Run44.
% t632_5: 63.2%% response time of the exact linearized channel.
% Differences are percentages relative to the linearized reference.
% No smoothing or new identification is applied to the noisy results.

gainDifferencePct = 100*(identified.K-Klinear)./Klinear;
timeDifferencePct = 100*(identified.tau-t632Linear)./t632Linear;

for scenario = 1:numel(scenarioNames)
    for k = 1:numel(stepSizes)
        fprintf('\n%s, +%.0f%% step\n',scenarioNames{scenario},100*stepSizes(k));
        fprintf('Channel       K4       K5    dK[%%]   tau4[s]  t632_5[s]   dt[%%]\n');
        for i = 1:2
            for j = 1:2
                fprintf('F%d->h%d  %8.5f %8.5f %8.2f %9.2f %10.2f %8.2f\n', ...
                    j,i,identified.K(i,j,k,scenario),Klinear(i,j), ...
                    gainDifferencePct(i,j,k,scenario), ...
                    identified.tau(i,j,k,scenario),t632Linear(i,j), ...
                    timeDifferencePct(i,j,k,scenario));
            end
        end
    end
end
