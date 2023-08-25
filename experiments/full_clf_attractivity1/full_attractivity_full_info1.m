%%  full_attractivity_full_info1.m
%   Description:
%       This script visualizes the trajectories that are meant to correctly slide our object in the world.
%       In order to safely and correctly slide, our system will need toreach the target set
%       (with velocity 0) for the full system's dynamics. In order to do this, it must first reach the
%       intermediate set where it will be released and will then reach the target.

clear all; close all; clc;

%% Add Functions to Path
addpath(genpath('../../src/'))

%% Constants
Theta_lb = 0.3; Theta_ub = 0.5;
Theta1 = Polyhedron('lb',Theta_lb,'ub',Theta_ub);
scalar_system = FrictionfulSliding('Theta', Theta1);

show_cylinder = false;

T_stop = 1;

X_Target_lb = 1.0; X_Target_ub = 1.3;
X_Target = Polyhedron('lb',X_Target_lb,'ub',X_Target_ub);
dtheta = 0.0125;

fs = 12;
x_lims = [0.0, 1.5];
xdot_lims = [1.0,2.0];
temp_axis = [x_lims, xdot_lims];

%% Plotting Cylinder
if show_cylinder
    [x1, y1, z1] = cylinder([0.1,0.5,0.8])
    figure;
    surf(x1,y1,z1)
end

%% Plotting Intersecting Cylinders
x_ub = [X_Target_ub, 0.0];
x_lb = [X_Target_lb, 0.0];
mu_k = 0.4;
g = 9.81;

dx = 0.001;
x = [x_lims(1):dx:x_lims(2)];
v = [xdot_lims(1):dx:xdot_lims(2)];
[x1pts, x2pts] = meshgrid(x, v);
x2 = []; y2 = []; z2 = [];

f_x = (1/(2*mu_k*g)) .* x2pts.^2 + x1pts;

in_range_fxs = (f_x >= x_lb(1)) & (f_x <= x_ub(1));
x2 = x1pts(in_range_fxs);
y2 = x2pts(in_range_fxs);

c_Xint = [0.9; 1.5];
Q = eye(2);

figure; %Plot the 2d slice
hold on;
scatter3(x2, y2, zeros(size(x2)) )
fcontour( ...
    @(x,xdot) ([x;xdot] - c_Xint )' * Q * ([x;xdot] - c_Xint), ...
    temp_axis, ...
    'LevelList',[0.01,0.025,0.05,0.1,0.15,0.2])

% axis(temp_axis)
view(2)

xlabel('$$x$$','Interpreter','latex','FontSize',fs)
ylabel('$$\dot{x}$$','Interpreter','latex','FontSize',fs)
% title(['Sampled Points That Reach Target set for $\mu_k=' num2str(mu_k) '$'],'Interpreter','latex')
% saveas(gcf,["X_int_slice_mu_0_4.png" ])