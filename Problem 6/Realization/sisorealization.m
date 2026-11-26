%
% sisorealization
%
%   Demonstrates state space realization of SISO 
%   transfer functions and computation of the
%   impulse and step responses.
%
%
clear
clc
close all

%            2
%  g(s) = ------- exp(-4*s)
%          s + 1
num = 2;
den = [1 1];
lambda = 4;


Ts = 0.01;              % Sampling time (or plot interval)            
tfinal = 10;            % Final plotting time
N = floor(tfinal/Ts);   % Last sample point computed

% Computes the impulse response
[h,th] = sisoctf2dimpulse(num,den,lambda,Ts,N);

% Computes the step response
[s,ts] = sisoctf2dstep(num,den,lambda,Ts,N);

% Plot the impulse response and the step response
figure(1)
subplot(211); plot(th,h); ylabel('Impulse response')
subplot(212); plot(ts,s); ylabel('Step response')
xlabel('time (min)')


%
% State space realization of transfer function with delay
% 
[a,b,c,d,tau] = sisoctf2css(num,den,lambda)
[ad,bd,cd,dd,l] = sisocss2dss(a,b,c,d,Ts,tau)

%
% Realization as a standard discrete time state space model
%
tol = 1.0e-12;
Nmax = 2*(lambda/Ts + length(den)-1);
[Ad,Bd,Cd,Dd]=mimoctf2dss({num},{den},lambda,Ts,Nmax,tol);
size(Ad,2)  

% -------------------------------------------------------------


%             (s-1)
% g1(s) = ------------ exp(-4*s)
%         (s+1)*(s+3) 
num1 = [1 -1];
den1 = conv([1 1],[1 3]);
lambda1 = 4;

Ts1 = 0.01;
tfinal1 = 10;
N1 = floor(tfinal1/Ts1);

[h1,th1] = sisoctf2dimpulse(num1,den1,lambda1,Ts1,N1);
[s1,ts1] = sisoctf2dstep(num1,den1,lambda1,Ts1,N1);

figure(2)
subplot(211); plot(th1,h1); ylabel('Impulse response')
subplot(212); plot(ts1,s1); ylabel('Step response')
xlabel('time (min)')


% -----------------------------------------------------------------

%
% g2(s) = exp(-4*s)
%
num2 = 1;
den2 = 1;
lambda2 = 4;

Ts2 = 0.01;
N2 = floor(10/Ts2);

[h2,th2] = sisoctf2dimpulse(num2,den2,lambda2,Ts2,N2);
[s2,ts2] = sisoctf2dstep(num2,den2,lambda2,Ts2,N2);

figure(3)
subplot(211); plot(th2,h2); ylabel('Impulse response')
subplot(212); plot(ts2,s2); ylabel('Step response')
xlabel('time (min)')


% ---------------------------------------------------------------

%
% g3(s) = 1
%
num3 = 1;
den3 = 1;
lambda3 = 0;

Ts3 = 0.01;
N3 = floor(10/Ts3);

[h3,th3] = sisoctf2dimpulse(num3,den3,lambda3,Ts3,N3);
[s3,ts3] = sisoctf2dstep(num3,den3,lambda3,Ts3,N3);

figure(4)
subplot(211); plot(th3,h3); ylabel('Impulse response')
subplot(212); plot(ts3,s3); ylabel('Step response')
xlabel('time (min)')
