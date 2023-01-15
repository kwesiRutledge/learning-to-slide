%%  example2_taylor
%   Description:
%       This script visualizes attempts to find an adaptive CBF for the
%       adaptive cruise control system that is discussed in the Taylor and
%       Ames paper "Adaptive Safety with Control Barrier Functions".

clear all; close all; clc;

%% Add Functions to Path
addpath(genpath('../../src/'))

%% Constants
acc1 = AdaptiveCruiseControl();

%% Plot Simple Trajectories of the System
x0 = [15;26];
tspan1 = [0,4];
[t1,y1] = ode45( ...
    @(t,x) acc1.dynamics(x,0), ...
    tspan1, x0 ...
    );

[t2,y2] = ode45( ...
    @(t,x) acc1.dynamics(x,-norm(x,2)*0.2), ...
    tspan1, x0 ...
    );

min_x1 = min([y1(:,1);y2(:,1)]);
min_x2 = min([y1(:,2);y2(:,2)]);
max_x1 = max([y1(:,1);y2(:,1)]);
max_x2 = max([y1(:,2);y2(:,2)]);

matching_axis1 = [ tspan1, min_x1, max_x1];

figure;
subplot(1,2,1)
plot(t1,y1(:,1))
xlabel('Time (s)')
ylabel('x')
axis(matching_axis1)
title('Zero Input Trajectory')

subplot(1,2,2)
plot(t2,y2(:,1))
xlabel('Time (s)')
ylabel('x')
axis(matching_axis1)
title('Slight corrective (u = -0.2x) trajectory')

%% Use the CLF-based controller for this work with known theta
%  ===========================================================

% Constants
alpha = 100;
dtheta = [1.0,0.1,0.05];
Theta = Polyhedron('lb',-dtheta,'ub',dtheta) + acc1.theta;

Gamma = eye(3);

deltaV = 0.1;

K_D = 10;
D_star = 30;

K_v = -1;
v_star = acc1.v0;

% Create a dummy program from the tutorial for SOSTools
syms v1 D1 f1 f2 f3 real;

u_min = -9.8*acc1.m;

Program1 = sosprogram([v1;D1]);
monom1 = monomials([v1;D1],[1:4]);
[Program1,l1] = sospolyvar(Program1,monom1);
[Program1,l2] = sospolyvar(Program1,monom1);
ha = alpha.^2 - (D1 - 1.8*v1 - alpha).^2;
Va_symb = (D1 - 1.8*v1).^2;

dV_dx_symb = gradient(Va_symb,[v1;D1])';

% Simulate
N_sims = 2;
T_sim = 20;
dt = 0.01;

decay_rate = 0;

Gamma_in = Gamma;

X0 = Polyhedron('lb',[10,20],'ub',[20,40]);

x0 = sampleFromPolytope(X0,N_sims);
th0 = sampleFromPolytope(acc1.Theta, N_sims);

[x, th0, V_history, th_history, u_history] = simulate_capa1_w_decay( ...
    acc1, ...
    Va_symb, ...
    T_sim, X0 , dt, ...
    Gamma_in, ...
    N_sims, decay_rate, ...
    {'v1','D1'}, {'f1','f2','f3'});

figure;
for x_index = 1:2
    subplot(2,1,x_index);
    hold on;
    for sim_index = 1:N_sims
        plot([0:T_sim],x{sim_index}(x_index,:))
        xlabel('Time Index (k)')
        ylabel(['$$x_' num2str(x_index) '(k)$$'],'Interpreter','latex')
    end
end

figure;
hold on;
for sim_index = 1:N_sims
    plot([0:T_sim],V_history{sim_index}(:))
    xlabel('Time Index (k)')
    ylabel(['$$V(x(k))$$'],'Interpreter','latex')
end