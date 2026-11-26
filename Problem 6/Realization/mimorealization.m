%
% mimorealization
%
%   Demonstrates state space realization of
%   a MIMO transfer function
%
clear
clc
close all

%
% Specify the system. Each SISO input-output function
% has the form g(s) = K/(T1+1)^2 * exp(-tau*s)
%
p = 4;          % number of outputs
m = 2;          % number of inputs 

tau = zeros(p,m);
T1 = zeros(p,m);
K = zeros(p,m);

K = [-0.2 -0.2;6 6.5; 23 33;0 0.95]
T1 = [15 15;15 14;15 12;0 3]
tau = [0 0; 15 0; 15 0;0 2]

%
% Compute cell arrays containing the transfer functions
%
num = cell(p,m);
den = cell(p,m);

for i=1:m
    for j=1:p
        num(j,i) = {K(j,i)};
        den(j,i) = {conv([T1(j,i) 1],[T1(j,i) 1])};
    end
end

%
% Compute a state space realization
%
Ts = 1.0;       % sampling time
nmax = 100;     % number of impulse response matrices computed for the realization
tol = 1.0e-6;   % tolerance (Hankel singular value cut-off)

[Ad,Bd,Cd,Dd]=mimoctf2dss(num,den,tau,Ts,nmax,tol);

%
% Compute the impulse response of the system
%

N = 150;                                  % Number of impulse response matrices
H1 = mimoctf2dimpulse(num,den,tau,Ts,N);  % Exact impulse response
H2 = mimodss2dimpulse(Ad,Bd,Cd,Dd,N);     % Impulse response of approximate state space model

figure(1)                                 % Plot the results (blue = exact, red = approx.)
hh1 = zeros(N+1,1);
hh2 = zeros(N+1,1);
for i=1:m
    for j=1:p
        for k=1:N+1
            hh1(k) = H1(j,i,k);
            hh2(k) = H2(j,i,k);
        end
        subplot(p,m,(j-1)*m+i)
        hold on;
        plot(hh1,'b-')
        plot(hh2,'r-');
        hold off
    end
end


%
% Compute the step response of the system
%

N = 150;                                  % Number of impulse response matrices
S1 = mimoctf2dstep(num,den,tau,Ts,N);     % Exact step response
S2 = mimodss2dstep(Ad,Bd,Cd,Dd,N);        % Step response of approximate state space model

figure(2)                                 % Plot the results (blue = exact, red = approx.)
ss1 = zeros(N+1,1);
ss2 = zeros(N+1,1);
for i=1:m
    for j=1:p
        for k=1:N+1
            ss1(k) = S1(j,i,k);
            ss2(k) = S2(j,i,k);
        end
        subplot(p,m,(j-1)*m+i)
        hold on;
        plot(ss1,'b-')
        plot(ss2,'r-');
        hold off
    end
end
