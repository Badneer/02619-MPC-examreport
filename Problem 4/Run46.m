clear; clc; close all;

scriptFolder = fileparts(mfilename('fullpath'));
data = load(fullfile(scriptFolder,'IdentifiedModel44.mat'),'G');
G = data.G;

Ts = 4;                   % Sampling period [s]
nMarkov = 300;            % Last lag, corresponding to 20 min
t = (0:nMarkov)*Ts;

% H(output, input, lag + 1)
H = zeros(2,2,nMarkov+1);

figure('Name','Discrete-time impulse response coefficients');
tiledlayout(2,2);

for i = 1:2
    for j = 1:2
        Gd = c2d(ss(G(i,j)),Ts,'zoh');
        [Ad,Bd,Cd,Dd] = ssdata(Gd);

        H(i,j,1) = Dd;              % H_0
        state = Bd;                 % x_1 after the unit input pulse

        for k = 1:nMarkov
            H(i,j,k+1) = Cd*state;  % H_k
            state = Ad*state;      % zero input after the pulse
        end

        nexttile
        stem(t,squeeze(H(i,j,:)),'.');
        title(sprintf('F_%d to h_%d',j,i));
        xlabel('Time after pulse [s]');
        ylabel(sprintf('(H_k)_{%d%d}',i,j));
        grid on
    end
end

KfromH = sum(H,3);
Kmodel = dcgain(G);

disp('Gain from the sum of Markov parameters')
disp(KfromH)

disp('Steady-state gain of the identified model')
disp(Kmodel)

fprintf('Maximum absolute gain difference: %.3e\n', ...
    max(abs(KfromH(:)-Kmodel(:))));