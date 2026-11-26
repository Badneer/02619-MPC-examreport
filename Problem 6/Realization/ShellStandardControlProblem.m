%
% ShellStandardControlProblem:
% ----------------------------
%  This script file demonstrates realization of the shell standard control
%  problem as discrete-time state space model using the realization
%  algorithms in the linear model predictive control toolbox.
%
%
%
clear
clc
close all

%
% Definition of the system by specification of gain, time contant, and time
% delay for each transfer function
%

K = [...
    4.05 1.77 5.88 1.20 1.44;...
    5.39 5.72 6.90 1.52 1.83;...
    3.66 1.65 5.53 1.16 1.27;...
    5.92 2.54 8.10 1.73 1.79;...
    4.13 2.38 6.23 1.31 1.26;...
    4.06 4.18 6.53 1.19 1.17;...
    4.38 4.42 7.20 1.14 1.26]

T = [...
    50 60 50 45 40;...
    50 60 40 25 20;...
     9 30 40 11  6;...
    12 27 20  5 19;...
     8 19 10  2 22;...
    13 33  9 19 24;...
    33 44 19 27 32]

tau = [...
      27 28 27 27 27;...
      18 14 15 15 15;...
       2 20  2  0  0;...
      11 12  2  0  0;...
       5  7  2  0  0;...
       8  4  1  0  0;...
       20 22 0  0  0]

%
% Compute the transfer functions
%
[p,m]=size(K);
num = cell(p,m);
den = cell(p,m);

for i=1:m
    for j=1:p
        num(j,i) = {K(j,i)};
        den(j,i) = {[T(j,i) 1]};
    end
end


%
% Compute a discrete-time state space realization
%
Ts = 1.0;          % Sampling time
nmax = 100;        % Maximum impulse responses used in state space realization
tol = 1.0e-6;      % tolerance 
tf = 300;          % final time for step response simulation

tic
[Ad,Bd,Cd,Dd,sH]=mimoctf2dss(num,den,tau,Ts,nmax,tol);
toc

%
% Plot the Hankel Singular Values
%
figure(1)
plot(log10(sH))
title('Hankel Singular Values')
xlabel('State dimension, n')
ylabel('log10(Hankel Singular Value)')


%
% Simulate step response of the true system 
% and the realized (approximative) system    
%

N = floor(tf/Ts);
[S1,ts]=mimoctf2dstep(num,den,tau,Ts,N);    % true system (red)
n1=size(Ad,1);
S2=mimodss2dstep(Ad,Bd,Cd,Dd,N);            % approx. system (blue)

% Plot the computed step responses
figure(2)
ss1 = zeros(N+1,1);
ss2 = zeros(N+1,1);
for i=1:m
    j=1;
        for k=1:N+1
            ss1(k) = S1(j,i,k);
            ss2(k) = S2(j,i,k);
        end
        subplot(p,m,(j-1)*m+i)
        title(strcat('u[',int2str(i),']'))
        hold on;
        plot(ts,ss1,'r-')
        plot(ts,ss2,'b-');
        hold off
    
    for j=2:p
        for k=1:N+1
            ss1(k) = S1(j,i,k);
            ss2(k) = S2(j,i,k);
        end
        subplot(p,m,(j-1)*m+i)
        hold on;
        plot(ts,ss1,'r-');  % true system
        plot(ts,ss2,'b-');  % approximative system
        hold off
    end
end

i=1;
for j=1:p
    subplot(p,m,(j-1)*m+i);
    hold on;
    ylabel(strcat('y[',int2str(j),']'));
    hold off;
end

% State dimension of realized system
n = size(Ad,2)