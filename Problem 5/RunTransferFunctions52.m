%% Problem 5.2 - Continuous-time transfer functions

% Compute the state-space models from Problem 5.1 first.

scriptFolder = fileparts(mfilename('fullpath'));
run(fullfile(scriptFolder,'RunLinearization51.m'));

%% 5.2 - Continuous-time transfer functions (zero initial deviations)
% G(s) = C*(s*I-A)^(-1)*B + D
% Columns 1:2: pump flows -> lower-tank levels.
% Columns 3:4: disturbance flows (Det/PC), or innovations (SDE) -> levels.
% Sampled measurement noise is additive and is not part of these dynamics.

GDet = tf(sysDet);
GPC = tf(sysPC);
GSDE = tf(sysSDE);

fprintf('5.2: deterministic transfer functions(GDet)\n');
display(GDet)
fprintf('5.2: piecewise stochastic transfer functions (GPC)\n');
display(GPC)
fprintf('5.2: SDE pump and innovation transfer functions (GSDE)\n');
display(GSDE)

