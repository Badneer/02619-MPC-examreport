%% Problem 5.1 - Continuous-time linearization
% All state-space models use deviations from the stated operating points.
% Outputs are the measured/controlled lower-tank levels [h1;h2].

%% Physical parameters and nominal pump flows

a = 1.2272*ones(4,1);
Atank = 380.1327*ones(4,1);
gamma = [0.6;0.7];
g = 981;
rho = 1.00;
p = [a;Atank;gamma;g;rho];
us = [300;300];

% Preserve the operating conditions defined in the Problem 2 scripts.
dsDet = [0;0];
dsPC = [20;20];
ad = [1/60;1/60];
stdF = [5;5];
sigmaF = stdF.*sqrt(2*ad);
Rvv = diag([0.1^2;0.1^2]);

%% 5.1 - Deterministic and piecewise-constant disturbance models
% delta_xdot = Ac*delta_x + Bc*delta_u + Ec*delta_d
% delta_y = Cc*delta_x             (deterministic)
% delta_y_k = Cc*delta_x_k + v_k   (stochastic), cov(v_k) = Rvv
% delta_z = Cc*delta_x for both models.

[AcDet,xsDet,hsDet] = linearizeTank(us,dsDet,p);
[AcPC,xsPC,hsPC] = linearizeTank(us,dsPC,p);

Bc = rho*[gamma(1) 0; 0 gamma(2); 0 1-gamma(2); 1-gamma(1) 0];
Ec = rho*[0 0; 0 0; 1 0; 0 1];
Cc = [1/(rho*Atank(1)) 0 0 0; 0 1/(rho*Atank(2)) 0 0];
Dc = zeros(2,4);

sysDet = ss(AcDet,[Bc Ec],Cc,Dc, ...
    'InputName',{'dF1','dF2','dF3','dF4'}, ...
    'OutputName',{'dh1','dh2'});
sysPC = ss(AcPC,[Bc Ec],Cc,Dc, ...
    'InputName',{'dF1','dF2','dF3','dF4'}, ...
    'OutputName',{'dh1','dh2'});

%% 5.1 - Augmented SDE model with Ornstein-Uhlenbeck disturbances
% xi = [m1;m2;m3;m4;F3;F4], xi_s = [xsPC;dsPC]
% d(delta_xi) = (AcSDE*delta_xi + BcSDE*delta_u)*dt + SigmaSDE*dW
% delta_y_k = CcSDE*delta_xi_k + v_k, cov(v_k) = Rvv
% delta_z = CcSDE*delta_xi

xis = [xsPC;dsPC];
AcSDE = [AcPC Ec; zeros(2,4) -diag(ad)];
BcSDE = [Bc;zeros(2,2)];
CcSDE = [Cc zeros(2,2)];
SigmaSDE = [zeros(4,2);diag(sigmaF)];

% The last two inputs represent the white-noise innovation channels.
% This ss object describes their shaping filter; simulation uses SigmaSDE*dW.
sysSDE = ss(AcSDE,[BcSDE SigmaSDE],CcSDE,Dc, ...
    'InputName',{'dF1','dF2','w3','w4'}, ...
    'OutputName',{'dh1','dh2'});

fprintf('5.1: nominal pump flows [cm^3/s]\n');
disp(us)
fprintf('Deterministic nominal disturbances [cm^3/s] and levels [cm]\n');
disp(dsDet)
disp(hsDet)
display(sysDet)
fprintf('Piecewise stochastic mean disturbances [cm^3/s] and levels [cm]\n');
disp(dsPC)
disp(hsPC)
display(sysPC)
fprintf('SDE model: augmented nominal state and continuous-time model\n');
disp(xis)
display(sysSDE)
fprintf('Stochastic sensor covariance Rvv [cm^2]\n');
disp(Rvv)

function [Ac,xs,hs] = linearizeTank(us,ds,p)
%linearizeTank Linearize the four mass balances at their nominal equilibrium.
    a = p(1:4);
    Atank = p(5:8);
    gamma = p(9:10);
    g = p(11);
    rho = p(12);

    q3s = (1-gamma(2))*us(2)+ds(1);
    q4s = (1-gamma(1))*us(1)+ds(2);
    qs = [gamma(1)*us(1)+q3s; gamma(2)*us(2)+q4s; q3s; q4s];
    hs = (qs./a).^2/(2*g);
    xs = rho*Atank.*hs;

    % Tank time constants
    Tau = (Atank./a).*sqrt(2*hs/g);
    % k_i = 1/Tau_i = a_i^2*g/(A_i*q_i,s).
    k = 1./Tau;
    Ac = [-k(1) 0 k(3) 0; 0 -k(2) 0 k(4); ...
          0 0 -k(3) 0; 0 0 0 -k(4)];
end
