function [s,ts] = sisodss2dstep(A,B,C,D,l,N,Ts)
%
% sisodss2dstep     Computes the discrete-time step response of a discrete-time state
%                   space system with delay.
%
%
%                   Discrete-time state space system with delay
%
%                       x[(k+1)*Ts] = A*x[k*Ts] + B*u[(k-l)*Ts]
%                          y[k*Ts]  = C*x[k*Ts] + D*u[(k-l)*Ts]
%
%                   The step response coefficients are
%
%                       s[k*Ts] = y[k*Ts]     k=0,1,...,N
%
%                   when 
%                       x[0] = 0  
%                       u[k*Ts] = 1     k = 0,1, ..., N-1
%
%
% Syntax: [s,ts] = sisodss2dstep(A,B,C,D,l,N,Ts)
%
%           A,B,C,D     :   Discrete-time state space matrices
%           l           :   Integer delay in discrete-time state space
%                           model
%           N           :   N*Ts is the final sampling time
%           Ts          :   Sampling time intervals
%
%           s           :   Step response coefficients
%           ts          :   Corresponding time vector
%


s = zeros(N+1,1);
ts = (0:Ts:N*Ts)';

[n,m]=size(B);
if n==0
        s(l+1:N+1) = D;    
else
    s(l+1) = D;
    x = B;
    for i=l+2:N+1
        s(i) = C*x+D;
        x = A*x+B;
    end
end
