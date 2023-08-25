%%  show_scalar_reachable_set
%   Description:
%       This script visualizes the set of states that can reach the target set
%       (with velocity 0) for the scalar system's dynamics.

clear all; close all; clc;

%% Add Functions to Path
addpath(genpath('../../src/'))

%% Constants
Theta_lb = -2.5; Theta_ub = -1.5;
Theta1 = Polyhedron('lb',Theta_lb,'ub',Theta_ub);
scalar_system = SimpleSystem1('Theta', Theta1);

T_stop = 1;

X_Target_lb = 0.2; X_Target_ub = 0.8;
X_Target = Polyhedron('lb',X_Target_lb,'ub',X_Target_ub);
dtheta = 0.0125;

x_dim = X_Target.Dim;
theta_dim = Theta1.Dim;

Gamma = 0.1;
dt = 0.02;

%% Compute PreSet
PreSet_XTheta = PolyUnion();
for theta = [Theta_lb: dtheta: Theta_ub]
    % Create a Polyhedron describing what initial conditions x0 will reach the
    % target set for this value of theta
    strip_i = Polyhedron( ...
        'lb',[X_Target_lb*exp(-T_stop*(1+theta)),        theta], ...
        'ub',[X_Target_ub*exp(-T_stop*(1+theta+dtheta)), theta+dtheta] ...
        );

    %Add strip to union
    PreSet_XTheta.add(strip_i);

end

%% Plotting PreSet with and without contours
temp_axis = [0.25,3.75,-2.6,-1.4];
fs = 22;

figure;
hold on;
plot(PreSet_XTheta, 'color','white') % Plot Pre Set
% plot([X_Target_lb,X_Target_lb],[Theta_lb,Theta_ub],':','Color','magenta','LineWidth',4.0)
% plot([X_Target_ub,X_Target_ub],[Theta_lb,Theta_ub],':','Color','magenta','LineWidth',4.0)
% plot(Polyhedron('lb',[X_Target_lb,Theta_lb],'ub',[X_Target_ub,Theta_ub]),'Color','magenta','alpha',0.5)

axis(temp_axis)

xlabel('$$x$$','Interpreter','latex','FontSize',fs)
ylabel('$$\theta$$','Interpreter','latex','FontSize',fs)

title(['Set of States that Reaches $$\mathcal{X}_T = [' num2str(X_Target_lb) ',' num2str(X_Target_ub) ']$$ at $$t=' num2str(T_stop) '$$'], ...
    'Interpreter','latex', ...
    'FontSize',fs ...
    )

scalar_R_fig_filename = ['original_target-T' num2str(T_stop) '.png'];
saveas(gcf,scalar_R_fig_filename)

[X_pts, Theta_pts] = meshgrid( ...
    [X_Target_lb-0.4:0.1:X_Target_ub+0.4], ...
    [Theta_lb:0.01:Theta_ub] ...
);

% Create a dummy program from the tutorial for SOSTools
syms x1 theta1 Va_cand(x1,theta1) real;
Va_cand = (x1 - (-2.5*theta1 - (13/4)) )^2

%% Compute Derivatives
state_names = {'x1'};
for x_index = 1:x_dim
    x_symb(x_index) = sym(state_names(x_index));
end
dVa_dx_symb = gradient(Va_cand,x_symb)';

theta_names = {'theta1'};
for theta_index = 1:theta_dim
    theta_symb(theta_index) = sym(theta_names(theta_index));
end
dVa_dth_symb = gradient(Va_cand,theta_symb)';

% Evaluate at each point
decay_rate = 0;
decrease_points = zeros(size(X_pts));
for row_index = [1:size(X_pts,1)]
    for col_index = [1:size(X_pts,2)]

        % Get x and theta
        x_k = X_pts(row_index, col_index);
        theta_hat_k = Theta_pts(row_index, col_index);

        % Evaluate barrier and clf at point
        Va_k = double( ...
                subs( ...
                    Va_cand, ...
                    state_parameter_pair_to_struct(x_k,theta_hat_k) ...
                ));
        dVa_k_dx = double( ...
            subs( ...
                dVa_dx_symb, ...
                state_parameter_pair_to_struct(x_k,theta_hat_k) ...
            ));

        % Evaluate the controller at each point
        [ u_k, optim_out ] = aCLF_control_capa2( ...
                x_k , ...
                Va_k, dVa_k_dx, dVa_dth_symb, ...
                scalar_system, ...
                decay_rate, Gamma, ...
                state_names, theta_names);

        % Compute next state
        x_kp1 = x_k + scalar_system.dynamics(x_k,u_k)*dt;
        theta_hat_kp1 = theta_hat_k + Gamma * aCLF_estimation_capa2(x_k,u_k,theta_hat_k,dVa_k_dx,scalar_system)*dt; 

        Va_kp1 = double( ...
                subs( ...
                    Va_cand, ...
                    state_parameter_pair_to_struct(x_kp1,theta_hat_kp1) ...
                ));

        if (Va_kp1 - Va_k)/dt < decay_rate * Va_k
            decrease_points(row_index, col_index) = 1.0;
        end

    end
end

figure;
hold on;
plot(PreSet_XTheta, 'color','white') % Plot Pre Set
surf(X_pts, Theta_pts, decrease_points-2.0)
xlabel('$$x$$','Interpreter','latex','FontSize',fs)
ylabel('$$\theta$$','Interpreter','latex','FontSize',fs)

axis(temp_axis)

level_set_fig_filename = ['target-set-with-level-set-lines.png'];
saveas(gcf,level_set_fig_filename)

%% Function for my Test aCLF for the scalar system

function [V_val] = aCLF_scalar1(x,theta_hat)
    V_val = (x - (-2.5*theta_hat - (13/4)))^2;
end
