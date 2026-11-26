function [Ad,Bd,Cd,Dd,sH]=mimodimpulse2dss(H,tol)
%
% mimodimpulse2dss  Computes a discrete-time state space model from
%                   discrete-time impulse response coefficients (Markov
%                   parameters) in the MIMO case.
%
%                   The discrete-time impulse response matrices of the MIMO
%                   system must be given as a (p-times-m-times-(N+1))
%                   matrix in which
%
%                       H(:,:,k)    is the impulse response matrix at time
%                                   k*Ts
%
%                   The discrete-time state space model is
%
%                       x[k+1] = Ad*x[k] + Bd*u[k]
%                         y[k] = Cd*x[k] + Dd*u[k]
%
%
%                   This model is realized by a balanced singular value
%                   decompositio of the Hankel matrix of the impulse
%                   response matrices. The number of impulse response
%                   matrices determines the computational load and also the
%                   maximum system order of the discrete time system that
%                   can be realized. Enough impulse response matrices must
%                   be included to determine the response of the system. If
%                   possible the number of impulse responses should be
%                   about 2 times the state dimension of the state space
%                   system.
%
%                   tol specifies the cut-off value of the singular values
%                   of the Hankel matrix of the impulse responses. These
%                   values, sH, are returned as well.
%
%                   The state space system realized is guaranteed to be
%                   both observable and controllable.
%
% Syntax: [Ad,Bd,Cd,Dd,sH]=mimodimpulse2dss(H,tol)
%
%               H           :   Sequence of impulse response matrices
%               tol         :   Hankel value tolerance cut-off
%
%               Ad,Bd,Cd,Dd :   Discrete-time system matrices
%               sH          :   Hankel matrix singular values
%

    


%
% Construct the Hankel Matrix
%
[p,m,Np1]=size(H);
N = floor((Np1-1)/2);

HH = zeros(p*(N+1),m*N);
for i=1:N
    kc1 = (i-1)*m+1;
    kc2 = kc1+m-1;
    for j=1:N+1
        kr1 = 1 + (j-1)*p;
        kr2 = kr1 + p - 1;
        HH(kr1:kr2,kc1:kc2) = H(:,:,2+(i-1)+(j-1));
    end
end

%
% SVD of the Hankel Matrix
%
[K,S,L] = svd(HH(1:p*N,1:m*N));


%
% Select the state dimension
%
sH = diag(S);
if nargin < 2
    tol = max(size(HH(1:p*N,1:m*N))') * max(sH) * eps;
end
nhat = sum(sH > tol);

%
% Compute the matrices
%
sqrtLambda = sqrt(S(1:nhat,1:nhat));
Bd = sqrtLambda*L(1:m,1:nhat)';
Cd = K(1:p,1:nhat)*sqrtLambda;
invsqrtLambda = diag(1./diag(sqrtLambda));
Ad = invsqrtLambda*K(:,1:nhat)'*HH(p+1:(N+1)*p,1:N*m)*L(:,1:nhat)*invsqrtLambda;

Dd = H(:,:,1);
