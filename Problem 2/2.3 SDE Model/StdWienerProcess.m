
function [W,T,dW] = StdWienerProcess(tf,N,nW,Ns,seed)
% STDWIENERPROCESS
% Generate independent standard Wiener processes

if nargin >= 5
    rng(seed);
end

dt = tf/N;

T = linspace(0,tf,N+1);

dW = sqrt(dt)*randn(nW,N,Ns);

W = cat(2,zeros(nW,1,Ns),cumsum(dW,2));

end
