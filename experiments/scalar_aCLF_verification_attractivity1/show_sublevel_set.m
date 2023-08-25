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

figure;
hold on;
plot(PreSet_XTheta, 'color','white') % Plot Pre Set
fcontour( ...
    @(x,theta_hat) (x - (-2.5*theta_hat - (13/4)))^2, ...
    temp_axis, ...
    'LevelList',[0.01,0.025,0.05,0.1,0.15,0.2])
xlabel('$$x$$','Interpreter','latex','FontSize',fs)
ylabel('$$\theta$$','Interpreter','latex','FontSize',fs)

axis(temp_axis)

level_set_fig_filename = ['target-set-with-level-set-lines.png'];
saveas(gcf,level_set_fig_filename)

%% Function for my Test aCLF for the scalar system

function [V_val] = aCLF_scalar1(x,theta_hat)
    V_val = (x - (-2.5*theta_hat - (13/4)))^2;
end
