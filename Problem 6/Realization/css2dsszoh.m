function [Phi,Gamma]=css2dsszoh(A,B,Ts)
%
% css2dsszoh Computes state space matrices for a zero-order hold discretization    
%
%           Let a continuous-time system be given by
%
%               d/dt x(t) = A*x(t) + B*u(t)
%
%           Assume that u(t) is implemented as a zero-order-hold
%
%               u(t) = u(k*Ts)          k*Ts <= t < (k+1)*Ts
%         
%           Then the discrete time equivalent system may be 
%           expressed as
%
%               x[(k+1)*Ts] = Phi*x[k*Ts] + Gamma*u[k*Ts]
%
%           The matrices Phi and Gamma are given by
%
%               Phi = exp(A*Ts)
%               Gamma = int_{0}^{Ts} exp(A*tau)*B dtau
%
%           Consequently, css2dsszoh may be regarded as
%           a function for computing these matrices/integrals
%
%
% Syntax:   [Phi,Gamma]=css2dsszoh(A,B,Ts)
%
%               A       :   Continuous-time matrix
%               B       :   Continuous-time matrix
%               Ts      :   Sampling time
%
%               Phi     :   Discrete-time matrix
%               Gamma   :   Discrete-time matrix
%

[n,m]=size(B);
M = [[A B]*Ts; zeros(m,n+m)];
S = expm(M);
Phi = S(1:n,1:n);
Gamma = S(1:n,n+1:n+m);

