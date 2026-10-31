
function xdot = FourTankSystem(t,x,u,d,p)
% FOURTANKSYSTEM
% Nonlinear dynamics of the Modified Four-Tank System
%
% dx/dt = f(t,x,u,d,p)
%
% States:
% x = [m1; m2; m3; m4]           [g]
%
% Manipulated inputs:
% u = [F1; F2]                   [cm^3/s]
%
% Disturbances:
% d = [F3; F4]                   [cm^3/s]
%
% Parameters:
% p = [a(1:4); A(1:4); gamma(1:2); g; rho]

%% Unpack parameters

a = p(1:4);
A = p(5:8);
gamma = p(9:10);
g = p(11);
rho = p(12);

%% Manipulated inputs and disturbances

F1 = u(1);
F2 = u(2);

F3 = d(1);
F4 = d(2);

%% Calculate liquid levels

h = x./(rho*A);

%% Outlet flows (Torricelli's law)

qout = a.*sqrt(2*g*h);

%% Inlet flows

qin = [
    gamma(1)*F1;
    gamma(2)*F2;
    (1-gamma(2))*F2 + F3;
    (1-gamma(1))*F1 + F4
];

%% Mass balances

xdot = zeros(4,1);

xdot(1) = rho*(qin(1) + qout(3) - qout(1));
xdot(2) = rho*(qin(2) + qout(4) - qout(2));
xdot(3) = rho*(qin(3) - qout(3));
xdot(4) = rho*(qin(4) - qout(4));

end
