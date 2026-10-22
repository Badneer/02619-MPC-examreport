function xdot = FourTankSystem(t,x,u,p)
% FOURTANKSYSTEM Model dx/dt = f(t,x,u,p) for 4-tank system

% Unpack states, MVs, and parameters
m = x;
F = u;
a = p(1:4,1);
A = p(5:8,1);
gamma = p(9:10,1);
g = p(11,1);
rho = p(12,1);

% Inflows
qin = zeros(4,1);
qin(1,1) = gamma(1)*F(1);
qin(2,1) = gamma(2)*F(2);
qin(3,1) = (1-gamma(2))*F(2);
qin(4,1) = (1-gamma(1))*F(1);

% Outflows
h = m./(rho*A);
qout = a.*sqrt(2*g*h);

% Differential equations
xdot = zeros(4,1);
xdot(1,1) = rho*(qin(1,1)+qout(3,1)-qout(1,1));
xdot(2,1) = rho*(qin(2,1)+qout(4,1)-qout(2,1));
xdot(3,1) = rho*(qin(3,1)-qout(3,1));
xdot(4,1) = rho*(qin(4,1)-qout(4,1));

end

