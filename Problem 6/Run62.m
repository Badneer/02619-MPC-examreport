%% Problem 6.2 - Dynamic vs steady-state KF, NO disturbance steps
% Reference truth: discretized augmented OU SDE from Problem 5.4.
clear;clc;close all;
scriptFolder=fileparts(mfilename('fullpath'));
dataFile=fullfile(scriptFolder,'Models61.mat');
assert(isfile(dataFile),'Run Run61.m first to create Models61.mat.');
load(dataFile,'AdSDE','BdSDE','CdSDE','QdSDE','Rdiscrete','Ts','xis');
A=AdSDE;B=BdSDE;C=CdSDE;Q=(QdSDE+QdSDE')/2;R=Rdiscrete;
nx=size(A,1);ny=size(C,1);
assert(nx==6 && ny==2,'Expected six augmented states and two measurements.');
N=round(20*60/Ts);t=(0:N)*Ts;
rng(100);


% ZOH pump deviations are zero. Unknown F3,F4 follow OU noise only.
u=zeros(2,N);


% Generate Gaussian increments even if Q is rank deficient.
[V,D]=eig(Q);L=V*diag(sqrt(max(diag(D),0)));
LR=chol(R,'lower');
X=zeros(nx,N+1);Y=zeros(ny,N+1);
for k=1:N
    X(:,k+1)=A*X(:,k)+B*u(:,k)+L*randn(nx,1);
end
for k=1:N+1
    Y(:,k)=C*X(:,k)+LR*randn(ny,1);
end
% Initial prior covariance: mass units [g], disturbance units [cm^3/s].
P0=diag([20^2*ones(1,4),5^2*ones(1,2)]);
x0=zeros(nx,1);

%% Dynamic (time-varying) KF
xPred=x0;Ppred=P0;
XhatDynamic=zeros(nx,N+1);
Kdynamic=zeros(nx,ny,N+1);
I=eye(nx);
for k=1:N+1
    S=C*Ppred*C'+R;
    K=(Ppred*C')/S;
    XhatDynamic(:,k)=xPred+K*(Y(:,k)-C*xPred);
    Kdynamic(:,:,k)=K;
    Ppost=(I-K*C)*Ppred*(I-K*C)'+K*R*K'; % Joseph form
    if k<=N
        xPred=A*XhatDynamic(:,k)+B*u(:,k);
        Ppred=A*Ppost*A'+Q;
    end
end

%% Static (steady-state-gain) KF
% Iterate the prior covariance Riccati recursion to convergence.
Pss=P0;
for iter=1:20000
    Sss=C*Pss*C'+R;
    Kss=(Pss*C')/Sss; % FILTERED gain, not one-step predictive gain
    Ppost=(I-Kss*C)*Pss*(I-Kss*C)'+Kss*R*Kss';
    Pnext=A*Ppost*A'+Q;
    if norm(Pnext-Pss,'fro')<1e-11*max(1,norm(Pss,'fro'))
        Pss=Pnext;
        break
    end
    Pss=Pnext;
end
assert(iter<20000,'Steady-state Riccati iteration did not converge.');
Kss=(Pss*C')/(C*Pss*C'+R);
XhatStatic=zeros(nx,N+1);xPred=x0;
for k=1:N+1
    XhatStatic(:,k)=xPred+Kss*(Y(:,k)-C*xPred);
    if k<=N
        xPred=A*XhatStatic(:,k)+B*u(:,k);
    end
end

%% Evaluate estimates: convert mass deviations into level deviations
Ztrue=C*X;Zdyn=C*XhatDynamic;Zstat=C*XhatStatic;
fprintf('Problem 6.2: no disturbance steps; Ts=%.1f s\n',Ts);

figure('Name','6.2 - Level estimation');
for j=1:2
    subplot(2,1,j);
    plot(t/60,Ztrue(j,:),'k','LineWidth',1.6);hold on
    plot(t/60,Y(j,:),'Color',[.7 .7 .7]);
    plot(t/60,Zdyn(j,:),'b','LineWidth',1.2);
    plot(t/60,Zstat(j,:),'r--','LineWidth',1.2);
    ylabel(sprintf('Delta h%d [cm]',j));grid on
    legend('True','Measured','Dynamic KF','Static KF','Location','best');
end
xlabel('Time [min]');

figure('Name','6.2 - Unknown disturbance estimation');
for j=1:2
    subplot(2,1,j);
    plot(t/60,xis(4+j)+X(4+j,:),'k','LineWidth',1.5);hold on
    plot(t/60,xis(4+j)+XhatDynamic(4+j,:),'b');
    plot(t/60,xis(4+j)+XhatStatic(4+j,:),'r--');
    ylabel(sprintf('F%d [cm^3/s]',j+2));grid on
    legend('True','Dynamic KF','Static KF','Location','best');
end
xlabel('Time [min]');

save(fullfile(scriptFolder,'Results62.mat'),'t','X','Y','XhatDynamic', ...
    'XhatStatic','Kdynamic','Kss','Ztrue','Zdyn','Zstat');
%% Problem 4 - Dynamic and static Kalman filters on the same measurements

identifiedData = load(dataFile,'sysIdentified_d');
assert(isfield(identifiedData,'sysIdentified_d'), ...
    'Run the updated Run61.m to save sysIdentified_d in Models61.mat.');

[Aid,Bid,Cid,Did] = ssdata(identifiedData.sysIdentified_d);

% The identified model has zero direct feedthrough.
assert(all(Did(:) == 0), ...
    'This filter block assumes zero direct feedthrough.');

nId = size(Aid,1);
Iid = eye(nId);

% Model the unmeasured disturbances as additive process noise.
% Map level-scale uncertainty into the identified state coordinates.
assert(rank(Cid) == size(Cid,1), ...
    'Expected two independent output rows in Cid.');
levelToState = pinv(Cid);

sigmaModel = 0.05;   % Assumed process increment scale [cm/sample]
sigmaInitial = 0.10; % Assumed initial level uncertainty [cm]

Qid = sigmaModel^2 * (levelToState*levelToState');
P0id = sigmaInitial^2 * (levelToState*levelToState');

% Steady-state prior covariance and filtering gain.
PssId = dare(Aid',Cid',Qid,R);
KssId = (PssId*Cid')/(Cid*PssId*Cid'+R);

% Both filters start at the nominal point in deviation coordinates.
xPredDynId = zeros(nId,1);
xPredStatId = zeros(nId,1);
PpredId = P0id;

ZdynId = zeros(size(Y));
ZstatId = zeros(size(Y));

for k = 1:size(Y,2)
    % Measurement update
    innovationCov = Cid*PpredId*Cid'+R;
    KdynId = (PpredId*Cid')/innovationCov;

    xDynId = xPredDynId + KdynId*(Y(:,k)-Cid*xPredDynId);
    xStatId = xPredStatId + KssId*(Y(:,k)-Cid*xPredStatId);

    ZdynId(:,k) = Cid*xDynId;
    ZstatId(:,k) = Cid*xStatId;

    % Time update
    if k <= size(u,2)
        xPredDynId = Aid*xDynId + Bid*u(:,k);
        xPredStatId = Aid*xStatId + Bid*u(:,k);

        correction = Iid-KdynId*Cid;
        PpostId = correction*PpredId*correction'+KdynId*R*KdynId';
        PpredId = Aid*PpostId*Aid'+Qid;
    end
end

%% Compare all four filters against the true levels

rmse = @(estimate) sqrt(mean((estimate-Ztrue).^2,2));

rmse62 = array2table( ...
    [rmse(Y), rmse(Zdyn), rmse(Zstat), rmse(ZdynId), rmse(ZstatId)], ...
    'VariableNames',{'Measured','P5_Dynamic','P5_Static', ...
    'P4_Dynamic','P4_Static'}, ...
    'RowNames',{'Tank1','Tank2'});

fprintf('\nProblem 6.2 - Level RMSE [cm]\n');
disp(rmse62);

figure('Name','6.2 - Level estimation with the Problem 4 model');
for j = 1:2
    subplot(2,1,j);
    plot(t/60,Y(j,:),'Color',[0.75 0.75 0.75]); hold on
    plot(t/60,Ztrue(j,:),'k','LineWidth',1.5);
    plot(t/60,ZdynId(j,:),'b','LineWidth',1.2);
    plot(t/60,ZstatId(j,:),'r--','LineWidth',1.2);
    ylabel(sprintf('Delta h%d [cm]',j));
    grid on
    legend('Measured','True','Dynamic KF (P4)','Static KF (P4)', ...
        'Location','best');
end
xlabel('Time [min]');

save(fullfile(scriptFolder,'Results62.mat'), ...
    'ZdynId','ZstatId','KssId','Qid','rmse62','-append');