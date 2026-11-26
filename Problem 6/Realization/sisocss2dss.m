function [Ad,Bd,Cd,Dd,l]=sisocss2dss(A,B,C,D,Ts,tau)
%
% sisocss2dss  Computes the discrete-time state space model of a
%              a continuous-time state space model with delay.
%
%              Continuous-time system
%
%                   d/dt x(t) = A*x(t) + B*u(t-tau)
%                        y(t) = C*x(t) + D*u(t-tau)
%
%               Zero-order-hold of the input
%
%                       u(t) = u(k*Ts)  k*Ts <= t < (k+1)*Ts
%
%               Discrete-time system
%
%                       x[(k+1)*Ts] = Ad*x[k*Ts] + Bd*u[(k-l)*Ts]
%                         y[k*Ts]   = Cd*x[k*Ts] + Dd*u[(k-l)*Ts]
%
%               The discrete time system matrices (Ad,Bd,Cd,Dd) are computed
%               assuming a zero-order-hold. l is the integer delay in the
%               discrete-time state space model
%
% Syntax: [Ad,Bd,Cd,Dd,l]=sisocss2dss(A,B,C,D,Ts,tau)
%
%           A,B,C,D     :   Continuous-time state space matrices
%           Ts          :   Sampling time
%           tau         :   Delay in continuous-time model
%
%           Ad,Bd,Cd,Dd :   Discrete-time state space matrices
%           l           :   Integer delay in discrete-time model
%



%
% Reference: Franklin, Powell & Workman: Digital Control of Dynamic
% Systems, 3rd edition, pp.110-114
%

if tau < 0
    error('Delay must be non-negative')    
end

[n,n]=size(A);
tmp = tau/Ts;
ll = ceil(tmp);
mm = ll - tmp;

if abs(mm) < eps
    % case: mm = 0
        if(n > 0)
            [Ad,Bd]=css2dsszoh(A,B,Ts);
            Cd = C; Dd = D;
            l = ll;
        else
        % n = 0
            Ad = []; Bd = [];
            Cd = []; Dd = D;
            l = ll;
        end
else
    % case: 0 < mm < 1
    if (n > 0)
        % Compute  Gamma1 = expm(A*m*Ts)*int_{0}{Ts-m*Ts} exp(A*eta) deta B
        %          Gamma2 = int_{0}{m*Ts} exp(A*eta) deta B
        tau2 = mm*Ts;
        tau1 = Ts-tau2;
        [N11,N12] = css2dsszoh(A,B,tau1);
        [M11,Gamma2] = css2dsszoh(A,B,tau2);
        Gamma1 = M11*N12;
        Phi = M11*N11;
        
        % Form the discrete time state space system
        Ad = [Phi Gamma1; zeros(1,n+1)];
        Bd = [Gamma2; 1];
        Cd = [C D];
        Dd = 0;
        l = ll-1;  % l = ll -1 = 1-1 = 0        
    else
        % n = 0
        Ad = [];
        Bd = [];
        Cd = [];
        Dd = D;
        l = ll;
    end    
end
    

