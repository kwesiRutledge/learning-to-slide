%% dai_p_toy_system1_functionalized.m
%Description:
%   Script meant to test the CLF verification function when you have a control affine system.
%   And input bounds.

clear all; close all; clc;

%% Add Functions to Path
addpath(genpath('../../src/'))

%% Constants

syms x1 x2 V(x1,x2) real

V_cand = x1^2 + x2^2;

ts = DaiP_ToySystem1();

%% Algorithm

sym_x = [x1;x2];

dim_x = length(sym_x);

[V, dV_dx, extraAlphas] = Verify_ca_CLF(sym_x, V_cand, ts)

%% Simulate
T_sim = 20;
N_sims = 10;
X0 = Polyhedron('lb',[-0.15,-0.15],'ub',[0.15,0.15]);
dt = 0.1;
[x,u_history, V_history] = simulate_ca_w_decay(ts, V, T_sim, X0 , dt, N_sims)

figure;
hold on;
for sim_index = 1:N_sims
    plot([0:T_sim],V_history{sim_index})
end
xlabel('k')
ylabel('$$V(x_k)$$','Interpreter','latex')

saveas(gcf,'daip_toysystem1_V_trajectories.png')

figure;
hold on;
for sim_index = 1:N_sims
    plot(x{sim_index}(1,:),x{sim_index}(2,:))
end
xlabel('$$(x_k)_1$$','Interpreter','latex')
ylabel('$$(x_k)_2$$','Interpreter','latex')

saveas(gcf,'daip_toysystem1_x_trajectories.png')

%% Function Definitions

function [ u_opt, optim_out ] = ca_CLF_output( x, V_x, dV_dx, ca_sys, decay_rate )
    %Description:
    %
    %Inputs
    %   V : Current value of the adaptive CLF
    %   dV_dx : Current value of the adaptive CLF's gradient with respect
    %           to the state x
    %           Should be a 1 x n vector where n is the dimension of the
    %           state
    %   


    % Constants
    settings_in = sdpsettings('verbose',0);
    dim_x = ca_sys.dim_x;
    dim_u = ca_sys.dim_u;

    if ~exist('decay_rate')
        decay_rate = 0;
    end

    decay_rate

    % Algorithm
    % =========

    u = sdpvar(dim_u,1)
    constraints = [];

    Lf_V = dV_dx * ca_sys.f(x);
    Lg_V = dV_dx * ca_sys.g(x);
    
    % Dynamics Constraint
    lyapunov_function_decrease_constraint = [Lf_V + Lg_V*u <= - decay_rate * V_x];
    
    constraints = constraints + lyapunov_function_decrease_constraint;

    % Input Bound Constraints
    input_constraint = [ca_sys.U.A * u <= ca_sys.U.b];
    constraints = constraints + input_constraint;

    % Optimize
    optim_out = optimize(constraints, u'*u ,settings_in);

    u_opt = value(u); % Save the feasible u

end

