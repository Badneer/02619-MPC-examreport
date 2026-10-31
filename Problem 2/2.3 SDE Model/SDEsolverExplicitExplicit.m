
function X = SDEsolverExplicitExplicit( ...
    ffun,gfun,T,x0,W,varargin)
% SDESOLVEREXPLICITEXPLICIT
% Euler-Maruyama method for stochastic systems

N = length(T);
nx = length(x0);

X = zeros(nx,N);
X(:,1) = x0;

for k = 1:N-1

    f = feval(ffun,T(k),X(:,k),varargin{:});
    sigma = feval(gfun,T(k),X(:,k),varargin{:});

    dt = T(k+1)-T(k);
    dW = W(:,k+1)-W(:,k);

    X(:,k+1) = X(:,k) + f*dt + sigma*dW;

end

end
