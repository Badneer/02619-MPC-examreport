function S = mimodss2dstep(Ad,Bd,Cd,Dd,N)
%
% mimodss2dstep     Computes the discrete-time step response matrices of
%                   a discrete-time MIMO state space system.
%
%                   Discrete-time system:
%
%                       x[k+1] = Ad*x[k] + Bd*u[k]
%                       y[k]   = Cd*x[k] + Dd*u[k]
%
%                   The step response matrices S(:,:,k) for this system are
%                   computed for k=0,1, ..., N.
%
%
% Syntax:   S = mimodss2dstep(Ad,Bd,Cd,Dd,N)
%
%           Ad,Bd,Cd,Dd :   Discrete-time state space matrices
%           N           :   Last sample to be included in the response
%
%           S           :   Step response matrices (p-times-m-times-(N+1))
%

H = mimodss2dimpulse(Ad,Bd,Cd,Dd,N);
[p,m,Nh]=size(H);
S = zeros(p,m,Nh);
S(:,:,1)=H(:,:,1);
for k=2:Nh
    S(:,:,k) = H(:,:,k)+S(:,:,k-1);
end
