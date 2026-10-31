
function sigma = FourTankDiffusion(t,x,u,p,ad,dbar,sigmaF)
% FOURTANKDIFFUSION
% Diffusion matrix for two independent
% Ornstein-Uhlenbeck disturbance processes.
%
% The stochastic noise enters F3 and F4 states.

sigma = zeros(6,2);

sigma(5,1) = sigmaF(1);
sigma(6,2) = sigmaF(2);

end
