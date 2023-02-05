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

fs = 22;

figure;
hold on;
plot(PreSet_XTheta, 'color','red') % Plot Pre Set
% plot([X_Target_lb,X_Target_lb],[Theta_lb,Theta_ub],':','Color','magenta','LineWidth',4.0)
% plot([X_Target_ub,X_Target_ub],[Theta_lb,Theta_ub],':','Color','magenta','LineWidth',4.0)
plot(Polyhedron('lb',[X_Target_lb,Theta_lb],'ub',[X_Target_ub,Theta_ub]),'Color','magenta','alpha',0.5)

xlabel('$$x$$','Interpreter','latex','FontSize',fs)
ylabel('$$\theta$$','Interpreter','latex','FontSize',fs)

title(['Set of States that Reaches $$\mathcal{X}_T = [' num2str(X_Target_lb) ',' num2str(X_Target_ub) ']$$ at $$t=' num2str(T_stop) '$$'], ...
    'Interpreter','latex', ...
    'FontSize',fs ...
    )

scalar_R_fig_filename = ['reach-T' num2str(T_stop) '.png'];
saveas(gcf,scalar_R_fig_filename)

%% I should verify this (only for lower left and upper left corners)
