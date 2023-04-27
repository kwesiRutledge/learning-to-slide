%%  show_sliding_problem_reachable_set
%   Description:
%       This script visualizes the set of states that can reach the target set
%       (with velocity 0) for the scalar system's dynamics.

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
x = [-3:dx:3];
v = [1:dx:2];
[x1pts, x2pts] = meshgrid(x, v);
x2 = []; y2 = []; z2 = [];

f_x = (1/(2*mu_k*g)) .* x2pts.^2 + x1pts;

in_range_fxs = (f_x >= x_lb(1)) & (f_x <= x_ub(1));
x2 = x1pts(in_range_fxs);
y2 = x2pts(in_range_fxs);

figure; %Plot the 2d slice
scatter(x2, y2)
title(['Sampled Points That Reach Target set for $\mu_k=' num2str(mu_k) '$'],'Interpreter','latex')
saveas(gcf,["X_int_slice_mu_0_4.png" ])

figure;
mu_k_list = [0.25,0.35,0.45];
for mu_k_index = [1:length(mu_k_list)]
    % Plot Each of these
    subplot(1,length(mu_k_list),mu_k_index)

    mu_k = mu_k_list(mu_k_index);

    f_x = (1/(2*mu_k*g)) .* x2pts.^2 + x1pts;

    in_range_fxs = (f_x >= x_lb(1)) & (f_x <= x_ub(1));
    x3 = x1pts(in_range_fxs);
    y3 = x2pts(in_range_fxs);

    scatter(x3, y3)

    title(['Intermediate Set for $\mu_k=' num2str(mu_k) '$' ],'Interpreter','latex')
    xlabel('x')
    ylabel('$\dot{x}$','Interpreter','latex')

    axis([0.0,1.5,1,2])

end
saveas(gcf,["X_int_slices.png"])