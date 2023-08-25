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

%% Verify Polynomial aCLF is valid
dim = 1;
Gamma = eye(dim);

% Create a dummy program from the tutorial for SOSTools
syms x theta_h Va_cand(x,theta_h) real;
Va_cand = (x - (-2.5*theta_h - (13/4)) )^2

%% Plot slices of this function
temp_axis = [-5,5,-70,100]

figure;
subplot(2,1,1)
ezplot(@(x) Va_fn(x,-2.4))
axis(temp_axis)

subplot(2,1,2)
ezplot(@(x) Va_fn(x,-1.1))
axis(temp_axis)

saveas(gcf,'images/reshaping_V.png')

function [Va_out] = Va_fn(x,theta_hat)
    Va_out = (4-x)*(x+4)*(x - (-2.5*theta_hat - (13/4)))^2;
end
