%% test_fn1_scalar_aclf_.m
%Description:
%   Script meant to test the aCLF verification function.

clear all; close all; clc;

%% Add Functions to Path
addpath(genpath('../../src/'))

%% Constants
Theta_lb = -2.5; Theta_ub = -1.5;
Theta1 = Polyhedron('lb',Theta_lb,'ub',Theta_ub);
scalar_system = SimpleSystem1(Theta1);

T_stop = 1;

%% Verify Polynomial aCLF is valid
dim = 1;
Gamma = eye(dim);

% Create a dummy program from the tutorial for SOSTools
syms x theta_h Va_cand(x,theta_h) real;
Va_cand = (x - (-2.5*theta_h - (13/4)) )^2

[Va, dVa_dx, dVa_dth] = Verify_aCLF( x, theta_h , Va_cand , scalar_system , Gamma )