%%  plot_V_true_vs_others.m
%   Description:
%       This script visualizes the value of V (given the true value of V).

clear all; close all; clc;

%% Add Functions to Path
addpath(genpath('../../src/'))

%% Constants
Theta_lb = -2.5; Theta_ub = -1.5;
Theta1 = Polyhedron('lb',Theta_lb,'ub',Theta_ub);
scalar_system = SimpleSystem1('Theta',Theta1);

T_stop = 1;

random_seed = 27;
rng(random_seed)

%% Verify Polynomial aCLF is valid
dim = 1;
Gamma = 0.5*eye(dim);

% Create a dummy program from the tutorial for SOSTools
syms x theta_h Va_cand(x,theta_h) real;
Va_cand = (x - (-2.5*theta_h - (13/4)) )^2

% [Va, dVa_dx, dVa_dth] = Verify_aCLF( x, theta_h , Va_cand , scalar_system , Gamma )

%% Simulate the System and Plot The Controller's performance

T_sim = 30;
N_sims = 3;
decay_rate = 0.5;
X0 = Polyhedron('lb',1.0,'ub',2.0);
dt = 0.1;

% Plot the lyapunov function over time.
[x, th0, Va_history, V_history, th_history, th_star, u_history] = simulate_capa2_with_V( ...
    scalar_system, ...
    Va_cand, ...
    T_sim, X0 , dt, ...
    Gamma, ...
    N_sims, decay_rate, ...
    {'x'}, {'theta_h'});

figure;
subplot(2,1,1)
hold on;
for sim_index = 1:N_sims
    plot([0:dt:T_sim*dt],x{sim_index})
end
xlabel('Time')
ylabel('x')
title('Unknown Theta, Exponentially Convergent aCLF Control')

subplot(2,1,2)
hold on;
for sim_index = 1:N_sims
    plot([0:dt:T_sim*dt],th_history{sim_index})
end
ylabel('$$\hat{\theta}$$','Interpreter','latex')
saveas(gcf,'images/scalar-capa2-x-theta-simulation.png')

figure;
hold on;
for sim_index = 1:N_sims
    plot([0:dt:T_sim*dt],V_history{sim_index})
end
xlabel('Time')
ylabel('V_t')
title('Unknown Theta, Exponentially Convergent aCLF Control')
saveas(gcf,'images/scalar-capa2-Va-simulation.png')

figure;
for sim_index = 1:N_sims
    subplot(N_sims,1,sim_index)
    hold on;
    plot([0:dt:T_sim*dt],V_history{sim_index}, ...
        'color','blue')
    plot([0:dt:T_sim*dt],Va_history{sim_index}, ...
        ':','color','black')
    xlabel('Time')
    legend('$$V$$','$$V_a$$','Interpreter','latex')
    title(['Simulation #' num2str(sim_index) ', V Trajectories'])
end
saveas(gcf,'images/scalar-capa2-Va-with-V.png')

figure;
hold on;
for sim_index = 1:N_sims
    theta_err = [];
    for tau = 0:T_sim
        theta_err = [theta_err, norm(th_history{sim_index}(:,tau+1) - th_star(:,sim_index))];
    end
    plot( ...
        [0:dt:T_sim*dt], ...
        theta_err ...
    )
end
xlabel('Time')
ylabel('$$\|\theta - \hat{\theta}\|$$','Interpreter','latex')
title('Parameter Estimation Error Over Time')
saveas(gcf,'images/scalar-capa2-theta-err.png')

%% Save data
date_string = datestr(now,'ddmmmyyyy-HHMM');
save(['data/scalar-capa2-full-aCLF-' date_string '.mat' ])
