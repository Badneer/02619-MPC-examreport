
function y = FourTankSensor(x,p)
% FOURTANKSENSOR
% Noise-free sensor function
%
% y = g(x,p)
%
% States:
% x = [m1; m2; m3; m4] [g]
%
% Measurements:
% y = [h1; h2]         [cm]
%
% Measurement noise is added separately:
% y_measured = g(x,p) + v

%% Parameters

A = p(5:8);
rho = p(12);

%% Tank levels

h = x./(rho*A);

%% Sensor outputs

y = h(1:2);

end
