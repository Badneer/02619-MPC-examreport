function [Ad,Bd,Cd,Dd,l]=sisoctf2dss(num,den,lambda,Ts)
%
% sisoctf2dss   Computes the discrete-time state space model with
%               integer delay for a continuous-time transfer function
%               with delay (SISO)
%
%               Let a SISO system be specified by the continuous-time transfer
%               function model
%
%                   y(s) = g(s)*u(s) 
%                        = h(s)*exp(-lambda*s)*u(s)
%                        = h(s)*u(s-lambda)
%
%                   h(s) = num(s)/den(s)
%
%                   den(s) = a(0)*s^na + a(1)*s^(na-1) + ... + a(na-1)*s + a(na)
%                   num(s) = b(0)*s^nb + b(1)*s^(nb-1) + ... + b(nb-1)*s + b(nb)
%
%               The time delay, lambda, must be non-negative and the transfer
%               function can only be realized if the denominator degree is
%               at least the same size as the numerator degree.
%
%               Zero-order-hold of the input (Ts: sampling time)
%
%                       u(t) = u(k*Ts)  k*Ts <= t < (k+1)*Ts
%
%               Discrete-time system
%
%                       x[(k+1)*Ts] = Ad*x[k*Ts] + Bd*u[(k-l)*Ts]
%                         y[k*Ts]   = Cd*x[k*Ts] + Dd*u[(k-l)*Ts]
%
%               The discrete time system matrices (Ad,Bd,Cd,Dd) are computed
%               assuming a zero-order-hold for the transfer function and sampling 
%               time Ts. l is the integer delay in the discrete-time state
%               space model
%
%
%
% Syntax: [Ad,Bd,Cd,Dd,l]=sisoctf2dss(num,den,lambda,Ts)
%
%           num         :   coefficients for the numerator in the transfer
%                           function
%           den         :   coefficients for the denominator in the transfer
%                           function
%           lambda      :   time delay (lambda >= 0)
%
%           Ad,Bd,Cd,Dd :   Discrete-time state space matrices
%           l           :   Integer delay in discrete-time model
%

%

[a,b,c,d,tau]=sisoctf2css(num,den,lambda);
[Ad,Bd,Cd,Dd,l]=sisocss2dss(a,b,c,d,Ts,tau);

