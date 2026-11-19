 %% Problem 6.3 - Kalman filtering with permanent disturbance steps

clear; clc; close all;

scriptFolder = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptFolder,'Models61.mat');

assert(isfile(dataFile), ...
    'Run the updated Run61.m first to create Models61.mat.');

data = load(dataFile,'sysPC_d','sysSDE_d','sysIdentified_d', ...
    'QdSDE','Rdiscrete','Ts','xis');

assert(isfield(data,'sysIdentified_d'), ...
    'Models61.mat must contain sysIdentified_d. Rerun the updated Run61.m.');

Ts = data.Ts;
R = data.Rdiscrete;
Q5 = (data.QdSDE + data.QdSDE')/2;

[A5,B5,C5,D5] = ssdata(data.sysSDE_d);
[A4,B4,C4,D4] = ssdata(data.sysIdentified_d);

assert(all(D5(:) == 0) && all(D4(:) == 0), ...
    'Expected zero feedthrough for both models.');
assert(data.sysSDE_d.Ts == Ts && data.sysPC_d.Ts == Ts ...
    && data.sysIdentified_d.Ts == Ts, ...
    'All models must use the same sample time.');

% sysPC_d inputs are [dF1 dF2 dF3 dF4].
% Columns 3 and 4 give the ZOH map for held disturbance flows.
% Add zero rows because the held offsets do not drive the OU fluctuations.
E5 = [data.sysPC_d.B(:,3:4); zeros(2,2)];

%% Generate one common linear test plant and measurement sequence

t = (0:round(20*60/Ts))*Ts;
N = numel(t)-1;
tStep = 8*60;
stepAmplitude = [10; 0];   % Permanent offsets in F3 and F4 [cm^3/s]

u = zeros(2,N);           % Nominal pump flows in deviation coordinates
flowOffset = stepAmplitude * double(t >= tStep);

rng(100);

[V,Lambda] = eig(Q5);
Lq = V*diag(sqrt(max(diag(Lambda),0)));
Lr = chol(R,'lower');

X = zeros(size(A5,1),N+1);

for k = 1:N
    X(:,k+1) = A5*X(:,k) + B5*u(:,k) ...
        + E5*flowOffset(:,k) + Lq*randn(size(A5,1),1);
end

Ztrue = C5*X;
Y = Ztrue + Lr*randn(2,N+1);

nominalFlows = data.xis(5:6);
nominalFlows = nominalFlows(:);
Ftrue = nominalFlows + X(5:6,:) + flowOffset;

%% Problem 5 - Integrating flow disturbances

n5 = size(A5,1);

A5aug = [A5 E5; zeros(2,n5) eye(2)];
B5aug = [B5; zeros(2,size(B5,2))];
C5aug = [C5 zeros(2,2)];

% Initial tuning for the flow-offset random walk [cm^3/s per increment].
sigmaFlowBias = 0.25;

Q5aug = blkdiag(Q5,sigmaFlowBias^2*eye(2));
P05 = diag([20^2*ones(1,4), 5^2*ones(1,2)]);
P05aug = blkdiag(P05,5^2*eye(2));

assert(rank([eye(n5+2)-A5aug; C5aug]) == n5+2, ...
    'Problem 5 augmented model failed the detectability check at z = 1.');

%% Problem 4 - Integrating output disturbances

n4 = size(A4,1);

A4aug = blkdiag(A4,eye(2));
B4aug = [B4; zeros(2,size(B4,2))];
C4aug = [C4 eye(2)];

% Reuse the Problem 4 process-noise assumption from Problem 6.2.
levelToState = pinv(C4);
Q4 = 0.05^2*(levelToState*levelToState');
P04 = 0.10^2*(levelToState*levelToState');

% Initial tuning for the output-offset random walk [cm per increment].
sigmaLevelBias = 0.02;

Q4aug = blkdiag(Q4,sigmaLevelBias^2*eye(2));
P04aug = blkdiag(P04,0.10^2*eye(2));

assert(rank([eye(n4+2)-A4aug; C4aug]) == n4+2, ...
    'Problem 4 augmented model failed the detectability check at z = 1.');

%% Run both filters for each augmented model

[Z5dyn,Z5stat,X5dyn,X5stat] = runKalmanPair( ...
    A5aug,B5aug,C5aug,Q5aug,R,P05aug,u,Y);

[Z4dyn,Z4stat] = runKalmanPair( ...
    A4aug,B4aug,C4aug,Q4aug,R,P04aug,u,Y);

% Baseline filters from the same model structures as in Problem 6.2.
% Rerun them using the new measurements with the permanent step.
[Z5baseDyn,Z5baseStat] = runKalmanPair(A5,B5,C5,Q5,R,P05,u,Y);
[Z4baseDyn,Z4baseStat] = runKalmanPair(A4,B4,C4,Q4,R,P04,u,Y);

%% Compare RMSE and mean estimation error

% Use the last five minutes to assess the remaining signed offset.
finalWindow = t >= t(end)-5*60;
assert(tStep < t(find(finalWindow,1)), ...
    'The final evaluation window must start after the disturbance step.');

estimates = {Y, Z5baseDyn, Z5baseStat, Z5dyn, Z5stat, ...
                Z4baseDyn, Z4baseStat, Z4dyn, Z4stat};

columnNames = {'Measured','P5_BaseDyn','P5_BaseStat', ...
    'P5_AugDyn','P5_AugStat','P4_BaseDyn','P4_BaseStat', ...
    'P4_AugDyn','P4_AugStat'};

rmse = zeros(2,numel(estimates));
meanError = zeros(2,numel(estimates));

for j = 1:numel(estimates)
    error = estimates{j}-Ztrue;
    rmse(:,j) = sqrt(mean(error.^2,2));
    meanError(:,j) = mean(error(:,finalWindow),2);
end

rmse63 = array2table(rmse, ...
    'VariableNames',columnNames,'RowNames',{'Tank1','Tank2'});

meanError63 = array2table(meanError, ...
    'VariableNames',columnNames,'RowNames',{'Tank1','Tank2'});

fprintf('\nProblem 6.3 - Level RMSE over the full run [cm]\n');
disp(rmse63);

fprintf('\nProblem 6.3 - Mean error over the last five minutes [cm]\n');
fprintf('Sign convention: estimated minus true level.\n');
disp(meanError63);

%% Plot level estimates and signed errors

plotLevelResults(t,tStep,'Problem 5',Ztrue,Y, ...
    Z5dyn,Z5stat,Z5baseDyn,Z5baseStat);

plotLevelResults(t,tStep,'Problem 4',Ztrue,Y, ...
    Z4dyn,Z4stat,Z4baseDyn,Z4baseStat);

% For P5, the total flow estimate includes the OU fluctuation and offset.
F5dyn = nominalFlows + X5dyn(5:6,:) + X5dyn(n5+(1:2),:);
F5stat = nominalFlows + X5stat(5:6,:) + X5stat(n5+(1:2),:);

figure('Name','6.3 - Unknown flow estimation with Problem 5');
for j = 1:2
    subplot(2,1,j);
    plot(t/60,Ftrue(j,:),'k','LineWidth',1.4); hold on
    plot(t/60,F5dyn(j,:),'b','LineWidth',1.2);
    plot(t/60,F5stat(j,:),'r--','LineWidth',1.2);
    xline(tStep/60,':','HandleVisibility','off');
    ylabel(sprintf('F%d [cm^3/s]',j+2));
    legend('True','Dynamic KF','Static KF','Location','best');
    grid on
end
xlabel('Time [min]');

save(fullfile(scriptFolder,'Results63.mat'), ...
    't','tStep','stepAmplitude','Ztrue','Y', ...
    'Z5dyn','Z5stat','Z4dyn','Z4stat', ...
    'Z5baseDyn','Z5baseStat','Z4baseDyn','Z4baseStat', ...
    'Ftrue','F5dyn','F5stat','rmse63','meanError63', ...
    'sigmaFlowBias','sigmaLevelBias');

%% Local functions

function [Zdyn,Zstat,Xdyn,Xstat] = runKalmanPair(A,B,C,Q,R,P0,u,Y)

    nx = size(A,1);
    numSamples = size(Y,2);
    I = eye(nx);

    % Solve the stationary Riccati equation for the prior covariance.
    Pss = dare(A',C',Q,R);
    Kss = (Pss*C')/(C*Pss*C'+R);

    xPredDyn = zeros(nx,1);
    xPredStat = zeros(nx,1);
    Ppred = P0;

    Xdyn = zeros(nx,numSamples);
    Xstat = zeros(nx,numSamples);

    for k = 1:numSamples
        Re = C*Ppred*C'+R;
        Kdyn = (Ppred*C')/Re;

        Xdyn(:,k) = xPredDyn + Kdyn*(Y(:,k)-C*xPredDyn);
        Xstat(:,k) = xPredStat + Kss*(Y(:,k)-C*xPredStat);

        if k < numSamples
            xPredDyn = A*Xdyn(:,k) + B*u(:,k);
            xPredStat = A*Xstat(:,k) + B*u(:,k);

            correction = I-Kdyn*C;
            Ppost = correction*Ppred*correction' + Kdyn*R*Kdyn';
            Ppred = A*Ppost*A'+Q;
            Ppred = (Ppred+Ppred')/2;
        end
    end

    Zdyn = C*Xdyn;
    Zstat = C*Xstat;
end

function plotLevelResults(t,tStep,modelName,Ztrue,Y, ...
                          Zdyn,Zstat,ZbaseDyn,ZbaseStat)

    figure('Name',['6.3 - ' modelName ' level estimates']);
    for j = 1:2
        subplot(2,1,j);
        plot(t/60,Y(j,:),'Color',[0.75 0.75 0.75]); hold on
        plot(t/60,Ztrue(j,:),'k','LineWidth',1.4);
        plot(t/60,Zdyn(j,:),'b','LineWidth',1.2);
        plot(t/60,Zstat(j,:),'r--','LineWidth',1.2);
        xline(tStep/60,':','HandleVisibility','off');
        ylabel(sprintf('Delta h%d [cm]',j));
        legend('Measured','True','Augmented dynamic','Augmented static', ...
            'Location','best');
        grid on
    end
    xlabel('Time [min]');

    figure('Name',['6.3 - ' modelName ' level estimation errors']);
    for j = 1:2
        subplot(2,1,j);
        plot(t/60,ZbaseDyn(j,:)-Ztrue(j,:),'Color',[0.35 0.35 0.35]); hold on
        plot(t/60,ZbaseStat(j,:)-Ztrue(j,:),'--','Color',[0.55 0.55 0.55]);
        plot(t/60,Zdyn(j,:)-Ztrue(j,:),'b','LineWidth',1.2);
        plot(t/60,Zstat(j,:)-Ztrue(j,:),'r--','LineWidth',1.2);
        xline(tStep/60,':','HandleVisibility','off');
        yline(0,':','HandleVisibility','off');
        ylabel(sprintf('Tank %d error [cm]',j));
        legend('Baseline dynamic','Baseline static', ...
            'Augmented dynamic','Augmented static','Location','best');
        grid on
    end
    xlabel('Time [min]');
end