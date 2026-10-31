
function z = FourTankOutput(x,p)
% FOURTANKOUTPUT
% Controlled output function
%
% z = h(x,p)
%
% States:
% x = [m1; m2; m3; m4] [g]
%
% Controlled outputs:
% z = [h1; h2]         [cm]

%% Parameters

A = p(5:8);
rho = p(12);

%% Tank levels

h = x./(rho*A);

%% Controlled outputs

z = h(1:2);

end
