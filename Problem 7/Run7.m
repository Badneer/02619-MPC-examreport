%% Problem 7 - Test the QP solver interface
clear; clc;

H = eye(2);
g = [-2; -2];

l = [0; 0];
u = [Inf; Inf];

A = [1 1];
bl = 1;
bu = 3;
xinit = [1; 1];

[x,info] = qpsolver(H,g,l,u,A,bl,bu,xinit);

disp('Optimal solution:')
disp(x)

disp('Solver exit flag:')
disp(info.exitflag)