function [s,ts]=sisoctf2dstep(num,den,lambda,Ts,N)
%
% sisoctf2dstep     The discrete time step response of a SISO
%                   continuous-time transfer function with delay
%
%                   Let a SISO system be specified by the continuous-time transfer
%                   function model
%
%                       y(s) = g(s)*u(s) 
%                            = h(s)*exp(-lambda*s)*u(s)
%                            = h(s)*u(s-lambda)
%
%                       h(s) = num(s)/den(s)
%
%                       den(s) = a(0)*s^na + a(1)*s^(na-1) + ... + a(na-1)*s + a(na)
%                       num(s) = b(0)*s^nb + b(1)*s^(nb-1) + ... + b(nb-1)*s + b(nb)
%
%                   The time delay, lambda, must be non-negative and the transfer
%                   function can only be realized if the denominator degree is
%                   at least the same size as the numerator degree.
%
%                   Zero-order-hold of the input (Ts: sampling time)
%
%                       u(t) = u(k*Ts)  k*Ts <= t < (k+1)*Ts
%                       
%                   The discrete-time step response coefficients
%
%                       s(k)    k = 0,1, ..., N
%
%                   are computed at times
%
%                       ts(k) = k*Ts
%
%                   assuming that u(t) is implemented as a zero-order-hold.
%
%
% Syntax: [s,ts]=sisoctf2dstep(num,den,lambda,Ts,N)
%
%           num         :   coefficients for the numerator in the transfer
%                           function
%           den         :   coefficients for the denominator in the transfer
%                           function
%           lambda      :   time delay (lambda >= 0)
%
%           s           :   Discrete time step response coefficients
%           ts          :   Corresponding time vector
%



[ad,bd,cd,dd,l]=sisoctf2dss(num,den,lambda,Ts);
[s,ts]=sisodss2dstep(ad,bd,cd,dd,l,N,Ts);
