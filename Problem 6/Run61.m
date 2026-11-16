%% Problem 6.1 - Discrete-time models from Problems 4 and 5

clear; clc; close all;

%% Load the discrete-time models from Problem 5.4

problem6Folder = fileparts(mfilename('fullpath'));
projectFolder = fileparts(problem6Folder);

discretizationScript = fullfile(projectFolder, ...
    'Problem 5','RunDiscretization54.m');

assert(isfile(discretizationScript), ...
    'Missing discretization script: %s',discretizationScript);

run(discretizationScript);

% Rebuild the paths because RunDiscretization54 runs in this workspace.
problem6Folder = fileparts(mfilename('fullpath'));
projectFolder = fileparts(problem6Folder);
problem4Folder = fullfile(projectFolder,'Problem 4');

%% Convert the identified model from Problem 4.4

identifiedFile = fullfile(problem4Folder,'IdentifiedModel44.mat');

assert(isfile(identifiedFile), ...
    'Missing %s. Run Problem 4/Run44.m first.',identifiedFile);

identified = load(identifiedFile,'G');

assert(isfield(identified,'G'), ...
    'IdentifiedModel44.mat must contain the transfer function matrix G.');

assert(identified.G.Ts == 0 && isequal(size(identified.G),[2 2]), ...
    'Expected G to be a continuous-time 2-by-2 transfer function matrix.');

sysIdentified_d = c2d(ss(identified.G),Ts,'zoh');

%% Report the discrete-time models

models = {sysDet_d,sysPC_d,sysSDE_d,sysIdentified_d};
names = {'Problem 5 - Deterministic', ...
         'Problem 5 - Piecewise stochastic', ...
         'Problem 5 - Augmented OU SDE', ...
         'Problem 4 - Identified from deterministic +10% steps'};

for k = 1:numel(models)
    model = models{k};
    fprintf('\n%s\n',names{k});
    fprintf('  %d states, %d inputs, %d outputs, Ts = %g s\n', ...
        size(model.A,1),size(model.B,2),size(model.C,1),model.Ts);
end

fprintf('\nProblem 4 - Discrete-time state-space realization:\n');
disp(sysIdentified_d);

fprintf('\nSDE process covariance Qd:\n');
disp(QdSDE);

fprintf('\nMeasurement covariance R:\n');
disp(Rdiscrete);

%% Save models for Problem 6.2

outputFile = fullfile(problem6Folder,'Models61.mat');

save(outputFile, ...
    'sysDet_d','sysPC_d','sysSDE_d','sysIdentified_d', ...
    'AdSDE','BdSDE','CdSDE','QdSDE','Rdiscrete','Ts','xis');

fprintf('\nSaved discrete-time models to:\n%s\n',outputFile);