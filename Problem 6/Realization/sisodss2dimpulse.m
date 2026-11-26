function [h,th] = sisodss2dimpulse(A,B,C,D,l,N,Ts)
%
% sisodss2dimpulse  Computes the discrete impulse response of a discrete-time state
%                   space system with delay (SISO).
%
%                   Discrete-time state space system with delay
%
%                       x[(k+1)*Ts] = A*x[k*Ts] + B*u[(k-l)*Ts]
%                          y[k*Ts]  = C*x[k*Ts] + D*u[(k-l)*Ts]
%
%                   The impulse response coefficients are
%
%                       h[k*Ts] = y[k*Ts]     k=0,1,...,N
%
%                   when 
%                       x[0] = 0 and 
%                       u[0] = 1 and 
%                       u[k*Ts] = 0 for k > 0
%
%
% Syntax: [h,th] = sisodss2dimpulse(A,B,C,D,l,N,Ts)
%
%           A,B,C,D     :   Discrete-time state space matrices
%           l           :   Integer delay in discrete-time state space
%                           model
%           N           :   N*Ts is the final sampling time
%           Ts          :   Sampling time intervals
%
%           h           :   Impulse response coefficients
%           th          :   Corresponding time vector
%

h = zeros(N+1,1);
th = (0:Ts:N*Ts)';

[n,m]=size(B);
if n==0
        h(l+1) = D;    
else
    h(l+1) = D;
    x = B;
    for i=l+2:N+1
        h(i) = C*x;
        x = A*x;
    end
end
