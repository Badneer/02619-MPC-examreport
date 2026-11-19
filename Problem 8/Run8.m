%% Problem 8 - Unconstrained MPC
scriptFolder = fileparts(mfilename('fullpath'));  % Problem 8
projectFolder = fileparts(scriptFolder);          % Repository root

modelFile = fullfile(projectFolder,'Problem 6','Models61.mat');

data = load(modelFile, ...
    'sysIdentified_d','AdSDE','BdSDE','CdSDE','Ts');
%% 1. Models and controller settings

[A4,B4,C4,D4] = ssdata(data.sysIdentified_d);
assert(all(D4(:) == 0), 'This prediction assumes D = 0.');

A5 = data.AdSDE;
B5 = data.BdSDE;
C5 = data.CdSDE;

Ts = data.Ts;
N = 60;                    % Prediction horizon in samples

Wz  = eye(2);               % Level tracking
Wu  = 0.01*eye(2);          % Deviation from preferred pump flows
Wdu = 0.10*eye(2);          % Flow changes between samples

mpc4 = mpcDesign(A4,B4,C4,N,Wz,Wu,Wdu);
mpc5 = mpcDesign(A5,B5,C5,N,Wz,Wu,Wdu);

%% 2. Test one regulator update at the nominal operating point

us = [300;300];             % Physical nominal pump flows [cm^3/s]
uPrevious = [0;0];          % Previous applied input in deviation units

r = [2;0];                  % Raise h1 by 2 cm and keep h2 at nominal
R = repmat(r,N,1);
Uref = zeros(2*N,1);        % Preferred pump flows are nominal

% Initial estimates for this test. In closed loop use the Kalman estimates.
xhat4 = zeros(size(A4,1),1);
xhat5 = zeros(size(A5,1),1);

[u4,U4,Z4,g4] = mpcCompute(mpc4,xhat4,uPrevious,R,Uref);
[u5,U5,Z5,g5] = mpcCompute(mpc5,xhat5,uPrevious,R,Uref);

fprintf('Prediction horizon: %d samples (%.0f s)\n',N,N*Ts);
disp('First physical pump commands [cm^3/s]')
disp(array2table([us+u4, us+u5], ...
    'VariableNames',{'Problem4','Problem5'}, ...
    'RowNames',{'F1','F2'}))

fprintf('Optimality residual, Problem 4: %.3e\n', ...
    norm(mpc4.H*U4+g4,inf));
fprintf('Optimality residual, Problem 5: %.3e\n', ...
    norm(mpc5.H*U5+g5,inf));

%% 3. Plot the predictions made at this update

tPrediction = (1:N)*Ts;
levels4 = reshape(Z4,2,N).';
levels5 = reshape(Z5,2,N).';

figure('Name','Problem 8 - MPC predictions');
tiledlayout(2,1);

for i = 1:2
    nexttile
    plot(tPrediction,levels4(:,i),'LineWidth',1.4);
    hold on
    plot(tPrediction,levels5(:,i),'--','LineWidth',1.4);
    yline(r(i),':','LineWidth',1.2);
    grid on
    xlabel('Time ahead [s]');
    ylabel(sprintf('Level deviation h_%d [cm]',i));
    legend('Problem 4','Problem 5','Setpoint','Location','best');
end

%% Local functions

function mpc = mpcDesign(A,B,C,N,Wz,Wu,Wdu)
% Build the prediction matrices and fixed terms of the regulator.

    nx = size(A,1);
    nu = size(B,2);
    nz = size(C,1);

    PhiX = zeros(N*nz,nx);
    Gamma = zeros(N*nz,N*nu);
    markov = zeros(nz,nu,N);
    Apower = eye(nx);

    for i = 1:N
        rows = (i-1)*nz + (1:nz);

        markov(:,:,i) = C*Apower*B;   % C*A^(i-1)*B
        Apower = Apower*A;
        PhiX(rows,:) = C*Apower;     % C*A^i
    end

    for i = 1:N
        rows = (i-1)*nz + (1:nz);

        for j = 1:i
            cols = (j-1)*nu + (1:nu);
            Gamma(rows,cols) = markov(:,:,i-j+1);
        end
    end

    % DeltaU = Lambda*U - I0*uPrevious.
    Lambda = kron(eye(N),eye(nu));
    if N > 1
        Lambda = Lambda - kron(diag(ones(N-1,1),-1),eye(nu));
    end
    I0 = [eye(nu); zeros((N-1)*nu,nu)];

    WzBar  = kron(eye(N),Wz);
    WuBar  = kron(eye(N),Wu);
    WduBar = kron(eye(N),Wdu);

    weightedGamma = WzBar*Gamma;
    weightedLambda = WduBar*Lambda;

    Hz  = weightedGamma.'*weightedGamma;
    Hu  = WuBar.'*WuBar;
    Hdu = weightedLambda.'*weightedLambda;

    mpc.H = Hz + Hu + Hdu;
    mpc.Mz = -weightedGamma.'*WzBar;
    mpc.Mu = -Hu;
    mpc.Mdu = -weightedLambda.'*WduBar*I0;

    mpc.PhiX = PhiX;
    mpc.Gamma = Gamma;
    mpc.nu = nu;
end

function [u,U,Z,g] = mpcCompute(mpc,xhat,uPrevious,R,Uref)
% Find the optimal input sequence and return its first move.

    [~,b] = mpcPredict(mpc,xhat,[]);
    c = R(:) - b;

    g = mpc.Mz*c + mpc.Mu*Uref(:) + mpc.Mdu*uPrevious(:);

    U = -(mpc.H\g);              % Solve H*U = -g
    u = U(1:mpc.nu);             % Apply only the first pair of commands
    Z = mpcPredict(mpc,xhat,U);
end

function [Z,b] = mpcPredict(mpc,xhat,U)
% Predict the stacked output deviations over the horizon.

    b = mpc.PhiX*xhat(:);

    Z = b;
    if ~isempty(U)
        Z = b + mpc.Gamma*U(:);
    end
end