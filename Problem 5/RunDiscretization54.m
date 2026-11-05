%% Problem 5.4 - Discrete-time state-space models

% ZOH discretization
% SDE process-noise covariance
% States, inputs and outputs are deviations from the operating points in 5.1.

scriptFolder = fileparts(mfilename('fullpath'));
run(fullfile(scriptFolder,'RunLinearization51.m'));

%% Sampling time
Ts = 4;  % [s], as in the course exercises and the Problem 2 scripts

%% Deterministic model
% delta_x(k+1) = AdDet*delta_x(k) + BdDet*delta_u(k) + EdDet*delta_d(k)
% delta_y(k) = CdDet*delta_x(k)
% ZOH: pump and disturbance flows are constant within each sample interval.

[AdDet,BallDet] = discretizeZOH(AcDet,[Bc Ec],Ts);
BdDet = BallDet(:,1:2);
EdDet = BallDet(:,3:4);
CdDet = Cc;
DdDet = Dc;
sysDet_d = ss(AdDet,[BdDet EdDet],CdDet,DdDet,Ts, ...
    'InputName',{'dF1','dF2','dF3','dF4'}, ...
    'OutputName',{'dh1','dh2'});

%% Piecewise-constant stochastic disturbance model
% delta_x(k+1) = AdPC*delta_x(k) + BdPC*delta_u(k) + EdPC*delta_d(k)
% delta_y(k) = CdPC*delta_x(k) + v(k), cov(v(k)) = Rvv
% delta_d(k) follows the disturbance process from Problem 2.2.
% A new independent disturbance is NOT assumed at every sample.

[AdPC,BallPC] = discretizeZOH(AcPC,[Bc Ec],Ts);
BdPC = BallPC(:,1:2);
EdPC = BallPC(:,3:4);
CdPC = Cc;
DdPC = Dc;
sysPC_d = ss(AdPC,[BdPC EdPC],CdPC,DdPC,Ts, ...
    'InputName',{'dF1','dF2','dF3','dF4'}, ...
    'OutputName',{'dh1','dh2'});

%% Augmented SDE model
% delta_xi(k+1) = AdSDE*delta_xi(k) + BdSDE*delta_u(k) + w(k)
% w(k) ~ N(0,QdSDE), independently over disjoint sample intervals.
% delta_y(k) = CdSDE*delta_xi(k) + v(k), cov(v(k)) = Rvv
% Only pump inputs are treated as ZOH inputs; diffusion is discretized below.

[AdSDE,BdSDE] = discretizeZOH(AcSDE,BcSDE,Ts);
CdSDE = CcSDE;
DdSDE = zeros(2,2);
sysSDE_d = ss(AdSDE,BdSDE,CdSDE,DdSDE,Ts, ...
    'InputName',{'dF1','dF2'},'OutputName',{'dh1','dh2'});

% QdSDE = integral_0^Ts exp(AcSDE*t)*SigmaSDE*SigmaSDE'*exp(AcSDE'*t) dt
nx = size(AcSDE,1);
Mnoise = [-AcSDE SigmaSDE*SigmaSDE'; zeros(nx,nx) AcSDE'];
PhiNoise = expm(Mnoise*Ts);
Phi12 = PhiNoise(1:nx,nx+1:2*nx);
Phi22 = PhiNoise(nx+1:2*nx,nx+1:2*nx);
QdSDE = Phi22'*Phi12;
QdSDE = (QdSDE+QdSDE')/2;  % Remove floating-point asymmetry.

% Sensor noise was already defined at sampling instants in Problem 2.
% Its covariance is unchanged; it is not multiplied by Ts.
Rdiscrete = Rvv;

%% Display the three discrete-time models
fprintf('\n5.4: discrete-time models, sampling time Ts = %.1f s\n',Ts);
display(sysDet_d)
display(sysPC_d)
display(sysSDE_d)
fprintf('SDE discrete process-noise covariance QdSDE:\n');
disp(QdSDE)
fprintf('Sampled measurement-noise covariance Rdiscrete (stochastic models):\n');
disp(Rdiscrete)

function [Ad,Bd] = discretizeZOH(Ac,Bc,Ts)
%discretizeZOH Exact ZOH discretization using the course block exponential.
    nx = size(Ac,1);
    nu = size(Bc,2);
    M = [Ac Bc; zeros(nu,nx) zeros(nu,nu)];
    Phi = expm(M*Ts);
    Ad = Phi(1:nx,1:nx);
    Bd = Phi(1:nx,nx+1:nx+nu);
end
