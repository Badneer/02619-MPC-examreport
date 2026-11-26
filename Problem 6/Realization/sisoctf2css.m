function [A,B,C,D,tau]=sisoctf2css(num,den,lambda)
%
% sisoctf2css   Realizes a continuous-time state space system
%               of a continuous-time transfer function (SISO)
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
%               The corresponding continuous-time state space realization
%               of this system is
%
%                  d/dt x(t) = A*x(t) + B*u(t-tau)
%                       y(t) = C*x(t) + D*u(t-tau)
%
%               sisoctf2css computes the matrices (A,B,C,D) and tau =
%               lambda.
%
%               The time delay must be non-negative and the transfer
%               function can only be realized if the denominator degree is
%               at least the same size as the numerator degree.
%
%               The model is realized in observer canonical form.
%
% Syntax: [A,B,C,D,tau]=sisoctf2css(num,den,lambda)
%
%           num         :   coefficients for the numerator in the transfer
%                           function
%           den         :   coefficients for the denominator in the transfer
%                           function
%           lambda      :   time delay (lambda >= 0)
%
%           (A,B,C,D)   :   state space matrices
%           tau         :   time delay (tau = lambda)
%       

[nden,nden2]=size(den);
[nnum,nnum2]=size(num);

if nargin==2
    lambda = 0;
end

if nden2 > 1
    nden = nden2;
    den = den';
end

if nnum2 > 1
    nnum = nnum2;
    num = num';
end


% Find first non-zero element in each transfer function
idden = min(find(abs(den)>eps));
den = den(idden:nden,1);
nden = length(den);

idnum=min( find( abs(num)>eps ) );
if isempty(idnum)
    num = 0;
    nnum = 1;
else
    num = num(idnum:nnum,1);
    nnum = length(num);
end

if abs(den(1,1)) < eps
    error('First element in denominator must be non-zero');
end

if nnum>nden
   error('Numerator degree must be less or equal to denominator degree');
elseif nnum<nden
    num = [zeros(nden-nnum,1); num];
end
    


n=nden-1;
if n > 0
    atilde = den(2:nden,1)/den(1,1);
    btilde0 = num(1,1)/den(1,1);
    btilde = num(2:nden,1)/den(1,1);
    
    % observer canonical form
    A = [-atilde(1:n-1,1) eye(n-1,n-1);-atilde(n,1) zeros(1,n-1)];
    B = btilde - atilde*btilde0;    
    C = zeros(1,n); C(1,1) = 1;
    D = btilde0;
    tau = lambda;
else
    A = [];
    B = [];
    C = [];
    D = num(1,1)/den(1,1);
    tau = lambda;
end






