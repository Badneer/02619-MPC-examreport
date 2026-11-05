%% Problem 5.5 - Markov parameters and comparison with Problem 4

% H0 = D, Hk = C*A^(k-1)*B for k >= 1.
% Each Hk is a 2x2 matrix: row = tank-level output, column = pump input.
% These are coefficients for a unit discrete pulse, with no division by Ts.

scriptFolder = fileparts(mfilename('fullpath'));
projectFolder = fileparts(scriptFolder);
problem4Folder = fullfile(projectFolder,'Problem 4');
identified = load(fullfile(problem4Folder,'IdentifiedModel44.mat'));
stepData = load(fullfile(problem4Folder,'StepData41.mat'));

run(fullfile(scriptFolder,'RunDiscretization54.m'));

% Run44 explicitly saves G for scenarioRef=1 (deterministic), kRef=1.
% Run41 uses us=[300;300], ds=[20;20] and the parameters listed below.

pProblem4 = [1.2272*ones(4,1);380.1327*ones(4,1);0.6;0.7;981;1.00];
if ~isequal(us,[300;300]) || ~isequal(dsPC,[20;20]) || ~isequal(p,pProblem4)
    error('The nominal conditions do not match the Problem 4 reference experiment.');
end
if any(abs(diff(stepData.T)-Ts) > 1e-10)
    error('Problems 4 and 5 must use the same sampling time for this comparison.');
end

% Use the recorded post-step horizon; do not invent an impulse horizon.
tMarkov = stepData.T(stepData.T >= stepData.tStep)-stepData.tStep;
nCoefficients = numel(tMarkov);

%% Pump-to-level Markov parameters for the three models from 5.4
HDet = markovParameters(AdDet,BdDet,CdDet,DdDet(:,1:2),nCoefficients);
HPC = markovParameters(AdPC,BdPC,CdPC,DdPC(:,1:2),nCoefficients);
HSDE = markovParameters(AdSDE,BdSDE,CdSDE,DdSDE,nCoefficients);

%% Discretize the actual identified model saved in Problem 4.4

% Compute its corresponding coefficients from the available G, using the
% same block-exponential ZOH method and sampling time as Problem 5.4.

sysIdentified_c = ss(identified.G);
A4 = sysIdentified_c.A;
B4 = sysIdentified_c.B;
C4 = sysIdentified_c.C;
D4 = sysIdentified_c.D;
nx4 = size(A4,1);
nu4 = size(B4,2);
Phi4 = expm([A4 B4; zeros(nu4,nx4) zeros(nu4,nu4)]*Ts);
Ad4 = Phi4(1:nx4,1:nx4);
Bd4 = Phi4(1:nx4,nx4+1:nx4+nu4);
sysIdentified_d = ss(Ad4,Bd4,C4,D4,Ts);
HIdentified = markovParameters(Ad4,Bd4,C4,D4,nCoefficients);

%% Display coefficients and matched operating-point comparison

fprintf('\n5.5: Markov parameters for pump inputs, Ts = %.1f s\n',Ts);
fprintf('H(1,:,:) stores H0; H(2,:,:) stores H1, etc.\n');
fprintf('Number of coefficients: %d, covering %.0f seconds.\n', ...
    nCoefficients,tMarkov(end));
fprintf('H0 = D = zero for all these strictly proper pump-to-level models.\n');
fprintf('Deterministic model H1 (nominal F3=F4=0):\n');
disp(squeeze(HDet(2,:,:)))
fprintf('Piecewise stochastic plant H1 (mean F3=F4=20):\n');
disp(squeeze(HPC(2,:,:)))
fprintf('Augmented SDE pump-channel H1 (mean F3=F4=20):\n');
disp(squeeze(HSDE(2,:,:)))
fprintf('Problem 4 identified reference H1 (deterministic +%.0f%%):\n', ...
    100*stepData.stepSizes(1));
disp(squeeze(HIdentified(2,:,:)))
fprintf('sysPC and sysSDE share the matched pump channels.\n');
fprintf('sysDet has a different operating point and is not used as a matched reference.\n');
fprintf('\nMaximum absolute coefficient differences from the identified model:\n');
maxMarkovDifference = zeros(2,2);
for i = 1:2
    for j = 1:2
        maxMarkovDifference(i,j) = max(abs(HPC(:,i,j)-HIdentified(:,i,j)));
        fprintf('F%d -> h%d: %.6g\n',j,i,maxMarkovDifference(i,j));
    end
end

%% Plot comparison of all four pump-to-level channels

figure('Name','Problem 5.5 - Markov parameter comparison');
for i = 1:2
    for j = 1:2
        subplot(2,2,(i-1)*2+j);
        plot(tMarkov/60,HPC(:,i,j),'b','LineWidth',1.5);
        hold on;
        plot(tMarkov/60,HIdentified(:,i,j),'r--','LineWidth',1.5);
        xlabel('Time after pulse [min]');
        ylabel('Markov coefficient [cm/(cm^3/s)]');
        title(sprintf('F%d to h%d',j,i));
        legend('Linearization (PC/SDE pump channels)', ...
            'Problem 4 identified reference','Location','best');
        grid on;
    end
end

function H = markovParameters(Ad,Bd,Cd,Dd,nCoefficients)
% markovParameters Compute the course discrete-time impulse coefficients.
    H = zeros(nCoefficients,size(Cd,1),size(Bd,2));
    H(1,:,:) = Dd;
    stateResponse = Bd;
    for k = 2:nCoefficients
        H(k,:,:) = Cd*stateResponse;
        stateResponse = Ad*stateResponse;
    end
end
