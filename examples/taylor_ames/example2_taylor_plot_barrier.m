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

%% Use the CLF-based controller for this work with known theta
%  ===========================================================

% Constants
alpha = 5;
dtheta = [1.0,0.1,0.05];
Theta = Polyhedron('lb',-dtheta,'ub',dtheta) + acc1.theta;

Gamma = eye(3);

%ha = alpha.^2 - (D1 - 1.8*v1 - alpha).^2;

figure;
ezsurf(@(v,D) ha_hybrid(v,D,alpha),[0,alpha,0,alpha])
colorbar()
hold on;
plot([-alpha:0.01:alpha],[-alpha:0.01:alpha]*(1.8)+0.1, ...
    'LineWidth',10.0,'Color','red')

view(0,90)

saveas(gcf,'images/barrier_plot.png')

function [ha_val] = ha_hybrid(v, D, alphaIn)
    if D - 1.8*v >= alphaIn
        ha_val = alphaIn.^2
    else
        ha_val = alphaIn^2 - (D - 1.8*v - alphaIn)^2
    end
end