function [x,info] = qpsolver(H,g,l,u,A,bl,bu,xinit)
%QPSOLVER Solve the convex quadratic program from Problem 7.
%
%   minimize    0.5*x'*H*x + g'*x
%   subject to  l <= x <= u
%               bl <= A*x <= bu
%
%   Use [] for an absent bound vector and +/-Inf for unbounded entries.

% 1. Prepare the objective and the initial guess.
g = g(:);
H = (H + H.')/2;
n = numel(g);

if isempty(xinit)
    xinit = zeros(n,1);
else
    xinit = xinit(:);
end

% 2. Convert two-sided linear bounds into Aineq*x <= bineq.
Aineq = [];
bineq = [];

if ~isempty(A)
    m = size(A,1);

    if isempty(bl)
        bl = -inf(m,1);
    end
    if isempty(bu)
        bu = inf(m,1);
    end

    bl = bl(:);
    bu = bu(:);

    hasUpper = bu < inf;
    hasLower = bl > -inf;

    Aineq = [ A(hasUpper,:);
        -A(hasLower,:)];

    bineq = [ bu(hasUpper);
        -bl(hasLower)];
end

% 3. Solve the QP. The active-set algorithm uses the initial guess.
options = optimoptions('quadprog', ...
    'Algorithm','active-set', ...
    'Display','off');

[x,fval,exitflag,output] = quadprog( ...
    H,g,Aineq,bineq,[],[],l(:),u(:),xinit,options);

% 4. Return diagnostics to the caller.
info.fval = fval;
info.exitflag = exitflag;
info.output = output;
end