
function f = FourTankDrift(t,x,u,p,ad,dbar,sigmaF)
% FOURTANKDRIFT
% Augmented 6-state modified four-tank system
%
% x(1:4) = Tank masses [g]
% x(5:6) = Disturbance flows F3,F4 [cm^3/s]
%
% ad     = Mean-reversion rates [1/s]
% dbar   = Mean disturbance flows [cm^3/s]
% sigmaF = Diffusion coefficients (used by diffusion function)

% Physical states
m = x(1:4);

% Current disturbance flows
d = x(5:6);

% Physical four-tank model from Common
fm = FourTankSystem(t,m,u,d,p);

% Ornstein-Uhlenbeck disturbance drift
fd = ad.*(dbar-d);

% Combined drift
f = [fm;fd];

end
