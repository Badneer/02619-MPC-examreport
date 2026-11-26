function H = mimodss2dimpulse(Ad,Bd,Cd,Dd,N)
%
% mimodss2dimpulse  Computes the discrete-time impulse response matrices of
%                   a discrete-time MIMO state space system.
%
%                   Discrete-time system:
%
%                       x[k+1] = Ad*x[k] + Bd*u[k]
%                       y[k]   = Cd*x[k] + Dd*u[k]
%
%                   The impulse response matrices H(:,:,k) for this system are
%                   computed for k=0,1, ..., N.
%
%
% Syntax:   H = mimodss2dimpulse(Ad,Bd,Cd,Dd,N)
%
%           Ad,Bd,Cd,Dd :   Discrete-time state space matrices
%           N           :   Last sample to be included in the response
%
%           H           :   Impulse response matrices (p-times-m-times-(N+1))
%

[p,m]=size(Dd);
H = zeros(p,m,N+1);

H(:,:,1) = Dd;
Gamma = Cd;
for k=2:N+1
    H(:,:,k) = Gamma*Bd;
    Gamma = Gamma*Ad;
end