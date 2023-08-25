%%scalar_lie_derivatives1_test.m
%Description
%   This script computes some of the Lie derivatives so that I can better
%   test some of these functions in the NeuralACLBF implementation.

%% Preamble
clear all; close all; clc;

%% Constants

% Get System
Theta_lb = 0.5; Theta_ub = 0.85;
Theta1 = Polyhedron('lb',Theta_lb,'ub',Theta_ub);
scalar_system = SimpleSystem1('Theta',Theta1);

x0 = 1.2;
h_theta_0 = scalar_system.theta;

dim_x = scalar_system.dim_x;
dim_u = scalar_system.dim_u;
dim_theta = scalar_system.dim_theta;

% Create aCLF
syms x1 theta1 a b c d;
Va_symb = x1.^2 + theta1.^2;

dVa_dx_symb  = gradient(Va_symb,[x1])';
dVa_dth_symb = gradient(Va_symb,[theta1])';

% Evaluate these at current x0, \hat{\theta}
Va = double( ...
    subs( ...
    Va_symb, state_parameter_pair_to_struct(x0,h_theta_0)...
    ));

dVa_dx = double( ...
    subs( ...
    dVa_dx_symb, state_parameter_pair_to_struct(x0,h_theta_0) ...
    ));

dVa_dth = double( ...
    subs( ...
        dVa_dth_symb, state_parameter_pair_to_struct(x0,h_theta_0) ...
    ));

%% Algorithm

Lf_V = dVa_dx * scalar_system.f(x0)
LF_V = dVa_dx * scalar_system.F(x0)
Lg_V = dVa_dx * scalar_system.g(x0)
LG_V = {};
G = scalar_system.G();
for theta_index = 1:dim_theta
    LG_V{theta_index} = dVa_dx * G{theta_index}
end

disp(['theta_hat = ' num2str(h_theta_0) ])
