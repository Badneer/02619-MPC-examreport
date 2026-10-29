
%% ---------------------------------------------------------
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

p = [a1;a2;a3;a4; A1;A2;A3;A4; ...
     gamma1;gamma2;g;rho];

%% ---------------------------------------------------------
% Disturbance variables
% ---------------------------------------------------------

d1 = 0;       % F3 [cm^3/s]
d2 = 0;       % F4 [cm^3/s]

d = [d1;d2];

%% ---------------------------------------------------------
% Nominal operating conditions
% ---------------------------------------------------------

F1_nom = 300;    % [cm^3/s]
F2_nom = 300;    % [cm^3/s]

us = [F1_nom;F2_nom];
ds = d;

%% ---------------------------------------------------------
% Simulation scenario
% ---------------------------------------------------------

t0 = 0;
tf = 20*60;

x0 = [0;0;0;0];

% Different F1 step sizes
step_percent = [0 5 10 25];

T_all = cell(length(step_percent),1);
X_all = cell(length(step_percent),1);

for k = 1:length(step_percent)

    F1 = F1_nom*(1 + step_percent(k)/100);

    u = [F1;F2_nom];

    % Simulate with disturbances
    [T_all{k},X_all{k}] = ode15s( ...
        @FourTankSystem,[t0 tf],x0,[],u,d,p);

end

%% ---------------------------------------------------------
% Compute nominal solution
% ---------------------------------------------------------

[T,X] = ode15s(@FourTankSystem, ...
               [t0 tf],x0,[],us,ds,p);

[nT,nX] = size(X);

a = p(1:4,1)';
A = p(5:8,1)';

% Compute liquid levels
H = zeros(nT,nX);

for i = 1:nT
    H(i,:) = X(i,:)./(rho*A);
end

% Compute outflows
Qout = zeros(nT,nX);

for i = 1:nT
    Qout(i,:) = a.*sqrt(2*g*H(i,:));
end

%% ---------------------------------------------------------
% Measurement and controlled output models
% ---------------------------------------------------------

% Measured variables: h1 and h2
Y = H(:,1:2);

% Controlled outputs: h1 and h2
Z = Y;

%% ---------------------------------------------------------
% Plot tank masses
% ---------------------------------------------------------

figure;

tank_order = [3 4 1 2];

for j = 1:4

    tank = tank_order(j);

    subplot(2,2,j)

    plot(T/60,X(:,tank)/1000, ...
         'b','LineWidth',1.5)

    xlabel('Time [min]')
    ylabel(sprintf('m_%d [kg]',tank))
    title(sprintf('Tank %d',tank))

    xlim([0 tf/60])
    grid on

end

%% ---------------------------------------------------------
% Compare different F1 step sizes
% ---------------------------------------------------------

figure;

for j = 1:4

    tank = tank_order(j);

    subplot(2,2,j)
    hold on

    for k = 1:length(step_percent)

        plot(T_all{k}/60, ...
             X_all{k}(:,tank)/1000, ...
             'LineWidth',1.5)

    end

    xlabel('Time [min]')
    ylabel(sprintf('m_%d [kg]',tank))
    title(sprintf('Tank %d',tank))

    legend('Nominal','+5%','+10%','+25%', ...
           'Location','best')

    xlim([0 tf/60])
    grid on
    hold off

end

%% ---------------------------------------------------------
% Linearization at steady state
% ---------------------------------------------------------

F1s = F1_nom;
F2s = F2_nom;

d1s = ds(1);
d2s = ds(2);

%% 1. Calculate steady-state flows

q3s = (1-gamma2)*F2s + d1s;
q4s = (1-gamma1)*F1s + d2s;

q1s = gamma1*F1s + q3s;
q2s = gamma2*F2s + q4s;

qs = [q1s;q2s;q3s;q4s];

%% 2. Calculate steady-state levels

a = [a1;a2;a3;a4];
A_tank = [A1;A2;A3;A4];

hs = (qs./a).^2/(2*g);

%% 3. Calculate steady-state masses

xs = rho*A_tank.*hs;

disp('Steady-state levels [cm]:')
disp(hs)

disp('Steady-state masses [kg]:')
disp(xs/1000)

%% ---------------------------------------------------------
% Linearized state-space matrices
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

%% B matrix - Manipulated variables

B_lin = rho*[gamma1       0;
                 0       gamma2;
                 0       1-gamma2;
              1-gamma1      0];

%% Bd matrix - Disturbance variables

B_dist = rho*[0 0;
              0 0;
              1 0;
              0 1];

%% C matrix

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

disp('B_dist = ')
disp(B_dist)

disp('C = ')
disp(C_lin)

disp('D = ')
disp(D_lin)

%% ---------------------------------------------------------
% Continuous-time state-space model
% ---------------------------------------------------------

sys_c = ss(A_lin,B_lin,C_lin,D_lin);

disp('Continuous-time state-space system:')
sys_c

%% ---------------------------------------------------------
% Discretize the linearized system
% ---------------------------------------------------------

Ts = 4;       % Sample time [s]

sys_d = c2d(sys_c,Ts,'zoh');

Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;
Dd = sys_d.D;

disp('Discrete-time A matrix:')
disp(Ad)

disp('Discrete-time B matrix:')
disp(Bd)

disp('Discrete-time C matrix:')
disp(Cd)

disp('Discrete-time D matrix:')
disp(Dd)

%% ---------------------------------------------------------
% Compare linear and nonlinear step responses
% ---------------------------------------------------------

t = (0:Ts:tf)';

stepPercent = [5 10 25];

figure

for j = 1:length(stepPercent)

    %% Step change in F1

    dF1 = F1_nom*stepPercent(j)/100;

    du = [dF1*ones(length(t),1), ...
          zeros(length(t),1)];

    %% Linear model

    dx0 = zeros(4,1);

    dy_lin = lsim(sys_d,du,t,dx0);

    hs_out = hs(1:2);

    y_lin = dy_lin + hs_out';

    %% Nonlinear model

    F1_step = F1_nom + dF1;

    u_nl = [F1_step;F2_nom];

    [Tnl,Xnl] = ode15s(@FourTankSystem, ...
                       [0 tf],xs,[],u_nl,ds,p);

    Hnl = Xnl./(rho*A_tank');

    %% Plot Tank 1

    subplot(2,3,j)

    plot(t/60,y_lin(:,1),'--','LineWidth',1.5)
    hold on

    plot(Tnl/60,Hnl(:,1),'-','LineWidth',1.5)

    xlabel('Time [min]')
    ylabel('h_1 [cm]')

    title(sprintf('Tank 1: +%d%% F_1',stepPercent(j)))

    legend('Linear','Nonlinear','Location','best')
    grid on

    %% Plot Tank 2

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
