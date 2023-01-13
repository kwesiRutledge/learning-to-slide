%% dai_p_toy_system1.m
%Description:
%   Script meant to test the CLF verification function when you have a .

clear all; close all; clc;

%% Add Functions to Path
addpath(genpath('../../src/'))

%% Constants

syms x1 x2 V(x1,x2) real

V_cand = x1^2 + x2^2;

ts = DaiP_ToySystem1();

%% Algorithm

V_U = ts.U.V;
n_VU = size(V_U,1);

sym_x = [x1;x2];

dim_x = length(sym_x);

rho0 = 0.1;

% Create Program
% ==============

syms offset0 real;
Program1 = sosprogram([sym_x]);

% Create Constraints
% ==================

% Create extra polynomial with adjustable parameter
offset_poly = -sym_x'*sym_x+offset0;

% Va should be p.d. in x, but not necessarily in theta

monom1 = monomials([sym_x],[0,4]);
[Program1, pd1] = sossosvar(Program1,monom1);

%Program1 = sosineq( Program1 , Va_cand - pd1 ); % Va_cand - pd1 >= 0

% Va satisfies gradient condition

for x_index = [1:dim_x]
    dV_dx_cand(x_index,1) = diff(V_cand,sym_x(x_index));
end

[Program1, lambda0] = sossosvar(Program1,monom1); %reate lambda 1
lambda_set = {};
for v_U_Index = [1:n_VU]
    [Program1, lambda_set{v_U_Index}] = sossosvar(Program1,monom1);
end

% Create gradient condition
grad_like_expression = (1+lambda0)*(V_cand - rho0)*(sym_x'*sym_x);
for v_index = [1:n_VU]
    % Get the theta at this vertex
    v_U = V_U(v_index,:)';

    grad_like_expression = grad_like_expression - ...
        lambda_set{v_index} * ( dV_dx_cand'*(ts.f(sym_x) + ts.g(sym_x)*v_U) );

end
% Enforce gradient is always decreasing.
Program1 = sosineq( Program1, grad_like_expression ) % grad_like_term >= 0 

% Solve
% =====

Program1 = sossolve( Program1 );
assert( (Program1.solinfo.info.pinf ~= 1) || (Program1.solinfo.info.dinf ~= 1) ); % Make sure that neither problem is infeasible?

V = sosgetsol(Program1,V_cand); %Getting solution for V
dV_dx = sosgetsol(Program1,dV_dx_cand); %Getting solution for dha_dx

%% Simulate
decay_rate = 1;

T_sim = 10; %Number of steps to simulate system forward.
dt = 0.1;
N_sims = 10;
x0 = unifrnd(-0.15,0.15,[2,N_sims]);

x = {}; u_history = {}; V_history = {};
for sim_index = 1:N_sims
    x_i = [x0(:,sim_index)];
    u_i = [];
    V_i = double(subs(V,struct('x1', x_i(1), 'x2', x_i(2))));

    for k = [0:T_sim-1]
        xi_k = x_i(:,k+1);
        
        Vi_k = double(subs(V,struct('x1', xi_k(1), 'x2', xi_k(2))));
        dVi_k_dx = double(subs(dV_dx,struct('x1', xi_k(1), 'x2', xi_k(2))))';

        [ui_k, optim_out] = ca_CLF_output(xi_k , Vi_k, dVi_k_dx, ts, decay_rate);
        xi_kp1 = xi_k + ts.dynamics(xi_k,ui_k)*dt;

        Vi_kp1 = double(subs(V,struct('x1', xi_kp1(1), 'x2', xi_kp1(2))));

        % Update the tracking
        x_i = [x_i, xi_kp1];
        u_i = [u_i, ui_k];
        V_i = [V_i, Vi_k];

    end

    x{sim_index} = x_i;
    u_history{sim_index} = u_i;
    V_history{sim_index} = V_i;

end

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

