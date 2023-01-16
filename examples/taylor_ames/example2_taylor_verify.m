%%  example2_taylor_VERIFY
%   Description:
%       This script visualizes attempts to find an adaptive CBF for the
%       adaptive cruise control system that is discussed in the Taylor and
%       Ames paper "Adaptive Safety with Control Barrier Functions".
%
%   Notes
%       Warning: The verification takes a long time for some reason (maybe
%       the polynomials are very large?)

clear all; close all; clc;

%% Add Functions to Path
addpath(genpath('../../src/'))

%% Constants
acc1 = AdaptiveCruiseControl();
acc2 = AdaptiveCruiseControlwNominal();

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
alpha = 5;
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
Va = (D1 - 1.8*v1).^2;
Va1 = D1.^2 + v1.^2;

Verify_aCLF_old([v1;D1],[f1;f2;f3],Va1, acc2, Gamma)
