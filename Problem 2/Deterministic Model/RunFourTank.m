% ---------------------------------------------------------
% Parameters
% ---------------------------------------------------------

a1 = 1.2272;
a2 = 1.2272;
a3 = 1.2272;
a4 = 1.2272;

A1 = 380.1327;
A2 = 380.1327;
A3 = 380.1327;
A4 = 380.1327;

gamma1 = 0.6;
gamma2 = 0.7;

g = 981;
rho = 1.00;

p = [a1;a2;a3;a4; A1;A2;A3;A4; gamma1;gamma2; g;rho];

% ---------------------------------------------------------
% Simulation scenario
% ---------------------------------------------------------

% ---------------------------------------------------------
% Simulation scenario
% ---------------------------------------------------------

t0 = 0;
tf = 20*60;

x0 = [0; 0; 0; 0];

F1_nom = 250;
F2 = 325;

% Different F1 step sizes
step_percent = [0 5 10 25];

% Store results
T_all = cell(length(step_percent),1);
X_all = cell(length(step_percent),1);

for k = 1:length(step_percent)

    % Change F1
    F1 = F1_nom*(1 + step_percent(k)/100);

    u = [F1; F2];

    % Simulate
    [T_all{k},X_all{k}] = ode15s( ...
        @FourTankSystem,[t0 tf],x0,[],u,p);

end

% ---------------------------------------------------------
% Compute the solution / Simulate
% ---------------------------------------------------------

[T,X] = ode15s(@FourTankSystem,[t0 tf],x0,[],u,p);

% Help variables
[nT,nX] = size(X);
a = p(1:4,1)';
A = p(5:8,1)';

% Compute liquid levels
H = zeros(nT,nX);

for i = 1:nT
    H(i,:) = X(i,:)./(rho*A);
end

% Compute flows out of each tank
Qout = zeros(nT,nX);

for i = 1:nT
    Qout(i,:) = a.*sqrt(2*g*H(i,:));
end

% --------------------------------------------------------
% Plot tank masses
% --------------------------------------------------------

figure;

% Tank 3
subplot(2,2,1)
plot(T/60, X(:,3), 'b', 'LineWidth', 1.5)
xlabel('time [min]')
ylabel('m_3 [kg]')
title('Tank 3')
xlim([0 tf/60])
grid off

% Tank 4
subplot(2,2,2)
plot(T/60, X(:,4), 'b', 'LineWidth', 1.5)
xlabel('time [min]')
ylabel('m_4 [kg]')
title('Tank 4')
xlim([0 tf/60])
grid off

% Tank 1
subplot(2,2,3)
plot(T/60, X(:,1), 'b', 'LineWidth', 1.5)
xlabel('time [min]')
ylabel('m_1 [kg]')
title('Tank 1')
xlim([0 tf/60])
grid off

% Tank 2
subplot(2,2,4)
plot(T/60, X(:,2), 'b', 'LineWidth', 1.5)
xlabel('time [min]')
ylabel('m_2 [kg]')
title('Tank 2')
xlim([0 tf/60])
grid off

figure

tank_order = [3 4 1 2];

for j = 1:4

    tank = tank_order(j);

    subplot(2,2,j)
    hold on

    for k = 1:length(step_percent)
        plot(T_all{k}/60, X_all{k}(:,tank), ...
            'LineWidth',1.5)
    end

    xlabel('Time [min]')
    ylabel(sprintf('m_%d [kg]',tank))
    title(sprintf('Tank %d',tank))

    legend('Nominal', '+5%', '+10%', '+25%', ...
        'Location','best')

    xlim([0 tf/60])
    grid on
    hold off

end

%% ---------------------------------------------------------
% Linearization of Four-Tank System at Steady State
% ---------------------------------------------------------

% Nominal inputs
F1s = 250;
F2s = 325;

us = [F1s; F2s];

%% 1. Calculate steady-state flows

q3s = (1-gamma2)*F2s;
q4s = (1-gamma1)*F1s;

q1s = gamma1*F1s + q3s;
q2s = gamma2*F2s + q4s;

qs = [q1s; q2s; q3s; q4s];


%% 2. Calculate steady-state levels

a = [a1; a2; a3; a4];
A_tank = [A1; A2; A3; A4];

hs = (qs./a).^2/(2*g);


%% 3. Calculate steady-state masses

xs = rho*A_tank.*hs;

disp('Steady-state levels [cm]:')
disp(hs)

disp('Steady-state masses [kg]:')
disp(xs)


%% ---------------------------------------------------------
% Linearized state-space matrices
% x = [m1 m2 m3 m4]'
% u = [F1 F2]'
% ---------------------------------------------------------

% Linearization constants
k1 = a1*g/(A1*q1s);
k2 = a2*g/(A2*q2s);
k3 = a3*g/(A3*q3s);
k4 = a4*g/(A4*q4s);


%% A matrix

A_lin = [ -k1     0      k3      0;
            0    -k2      0      k4;
            0     0     -k3      0;
            0     0       0     -k4];


%% B matrix

B_lin = rho*[ gamma1       0;
                 0       gamma2;
                 0       1-gamma2;
              1-gamma1      0];


%% C matrix
% Outputs are h1 and h2

C_lin = [1/(rho*A1)      0        0   0;
              0      1/(rho*A2)   0   0];


%% D matrix

D_lin = zeros(2,2);


%% ---------------------------------------------------------
% Display matrices
% ---------------------------------------------------------

disp('A = ')
disp(A_lin)

disp('B = ')
disp(B_lin)

disp('C = ')
disp(C_lin)

disp('D = ')
disp(D_lin)


%% ---------------------------------------------------------
% Create continuous-time state-space model
% ---------------------------------------------------------

sys_c = ss(A_lin,B_lin,C_lin,D_lin);

disp('Continuous-time state-space system:')
sys_c




%% ---------------------------------------------------------
% Discretize the linearized system
% ---------------------------------------------------------

Ts = 4;       % Sample time [s]

sys_d = c2d(sys_c, Ts, 'zoh');

% Extract discrete-time matrices
Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;
Dd = sys_d.D;

%% Display results

disp('Discrete-time A matrix:')
disp(Ad)

disp('Discrete-time B matrix:')
disp(Bd)

disp('Discrete-time C matrix:')
disp(Cd)

disp('Discrete-time D matrix:')
disp(Dd)

%% =========================================================
% Compare linear and nonlinear step responses
% ==========================================================

% Simulation time
tf = 20*60;                 % [s]
t = (0:Ts:tf)';             % use same Ts as discrete model

% Nominal inputs
F1s = 250;
F2s = 325;
us = [F1s; F2s];

% Step sizes
stepPercent = [5 10 25];

figure

for j = 1:length(stepPercent)

    %% -----------------------------------------------------
    % Step change in F1
    % ------------------------------------------------------

    dF1 = F1s*stepPercent(j)/100;

    % Input deviation for linear model
    % F1 changes, F2 remains unchanged
    du = [dF1*ones(length(t),1), ...
          zeros(length(t),1)];


    %% -----------------------------------------------------
    % LINEAR MODEL
    % ------------------------------------------------------

    % Linear model works in deviation variables
    % Initial deviation from steady state = 0

    dx0 = zeros(4,1);

    % Use discrete linear model
    dy_lin = lsim(sys_d,du,t,dx0);

    % Nominal steady-state outputs
    hs_out = [hs(1); hs(2)];

    % Convert deviation output to actual tank level
    y_lin = dy_lin + hs_out';


    %% -----------------------------------------------------
    % NONLINEAR MODEL
    % ------------------------------------------------------

    % Actual input after step
    F1_step = F1s + dF1;

    u_nl = [F1_step; F2s];

    % Start nonlinear model at nominal steady state
    [Tnl,Xnl] = ode15s(@FourTankSystem, ...
                       [0 tf], ...
                       xs, ...
                       [], ...
                       u_nl,p);

    % Convert nonlinear masses to levels
    Hnl = zeros(length(Tnl),4);

    for i = 1:length(Tnl)
        Hnl(i,:) = Xnl(i,:)./(rho*A_tank');
    end


    %% -----------------------------------------------------
    % Plot Tank 1
    % ------------------------------------------------------

    subplot(2,3,j)

    plot(t/60,y_lin(:,1),'--','LineWidth',1.5)
    hold on

    plot(Tnl/60,Hnl(:,1),'-','LineWidth',1.5)

    xlabel('Time [min]')
    ylabel('h_1 [cm]')

    title(sprintf('Tank 1: +%d%% F_1',stepPercent(j)))

    legend('Linear','Nonlinear','Location','best')

    grid on


    %% -----------------------------------------------------
    % Plot Tank 2
    % ------------------------------------------------------

    subplot(2,3,j+3)

    plot(t/60,y_lin(:,2),'--','LineWidth',1.5)
    hold on

    plot(Tnl/60,Hnl(:,2),'-','LineWidth',1.5)

    xlabel('Time [min]')
    ylabel('h_2 [cm]')

    title(sprintf('Tank 2: +%d%% F_1',stepPercent(j)))

    legend('Linear','Nonlinear','Location','best')

    grid on

end
%% ---------------------------------------------------------
% Continuous-time transfer functions
% ---------------------------------------------------------

clear tf

G = tf(sys_c);

G11 = G(1,1);    % F1 -> h1
G12 = G(1,2);    % F2 -> h1
G21 = G(2,1);    % F1 -> h2
G22 = G(2,2);    % F2 -> h2

disp('F1 -> h1')
G11

disp('F2 -> h1')
G12

disp('F1 -> h2')
G21

disp('F2 -> h2')
G22