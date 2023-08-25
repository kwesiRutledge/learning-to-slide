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
scalar_system = SimpleSystem1(Theta1);

T_stop = 1;

X_Target_lb = 0.2; X_Target_ub = 0.8;
X_Target = Polyhedron('lb',X_Target_lb,'ub',X_Target_ub);
dtheta = 0.0125;
PreSet_XTheta = PolyUnion();
PreSet_XTheta_3d = PolyUnion();
for theta = [Theta_lb: dtheta: Theta_ub]
    % Create a Polyhedron describing what initial conditions x0 will reach the
    % target set for this value of theta
    strip_i = Polyhedron( ...
        'lb',[X_Target_lb*exp(-T_stop*(1+theta)),        theta], ...
        'ub',[X_Target_ub*exp(-T_stop*(1+theta+dtheta)), theta+dtheta] ...
        );

    %Add strip to union
    PreSet_XTheta.add(strip_i);
    PreSet_XTheta_3d.add(strip_i * Polyhedron('lb',20,'ub',21));

end

fs = 22;

figure;
hold on;
plot(PreSet_XTheta) % Plot Pre Set
% plot([X_Target_lb,X_Target_lb],[Theta_lb,Theta_ub],':','Color','magenta','LineWidth',4.0)
% plot([X_Target_ub,X_Target_ub],[Theta_lb,Theta_ub],':','Color','magenta','LineWidth',4.0)
plot(Polyhedron('lb',[X_Target_lb,Theta_lb],'ub',[X_Target_ub,Theta_ub]),'Color','magenta','alpha',0.5)

xlabel('$$x$$','Interpreter','latex','FontSize',fs)
ylabel('$$\theta$$','Interpreter','latex','FontSize',fs)

title(['Set of States that Reaches $$\mathcal{X}_T = [' num2str(X_Target_lb) ',' num2str(X_Target_ub) ']$$ at $$t=' num2str(T_stop) '$$'], ...
    'Interpreter','latex', ...
    'FontSize',fs ...
    )

scalar_R_fig_filename = ['intermediate-T' num2str(T_stop) '.png'];
saveas(gcf,scalar_R_fig_filename)

%% Verify Polynomial aCLF is valid
dim = 1;
Gamma = eye(dim);

% Create a dummy program from the tutorial for SOSTools
syms x theta_h Va(x,theta_h) real;

Program1 = sosprogram([x;theta_h]);
monom1 = monomials([x;theta_h],[1:4]);
[Program1,l1] = sospolyvar(Program1,monom1);
[Program1,l2] = sospolyvar(Program1,monom1);
Va = (x - (-2.5*theta_h - (13/4)) )^2;

monom2 = monomials([x;theta_h],[0:2]);
[Program1,u1] = sospolyvar(Program1,monom2);

dVa_dx(1,1) = diff(Va,x);

dVa_dth(1,1) = diff(Va,theta_h);

% Create constraints
u_min = -1;
u_max =  1;
tempExpr =  dVa_dx'* ( scalar_system.f(x) + scalar_system.F(x)* ( Theta_lb - Gamma * dVa_dth ) + scalar_system.g(x)* u_min) - l1 ;
Program1 = sosineq( Program1 , tempExpr ); %F3

tempExpr =  dVa_dx'* ( scalar_system.f(x) + scalar_system.F(x)* ( Theta_ub - Gamma * dVa_dth ) + scalar_system.g(x)* u_min) - l2 ;
Program1 = sosineq( Program1 , tempExpr ); %F3

%Both helper functions are p.d.
Program1 = sosineq( Program1 , l1 ); 
Program1 = sosineq( Program1 , l2 );

% Optimize
Program1 = sossolve( Program1 );
assert( (Program1.solinfo.info.pinf ~= 1) || (Program1.solinfo.info.dinf ~= 1) ); % Make sure that neither problem is infeasible?

% Plot the obtained lyapunov function
SOL_Va = sosgetsol(Program1,Va); %Getting solution for V
SOL_dVa_dx = sosgetsol(Program1,dVa_dx); %Getting solution for dha_dx
SOL_dVa_dth = sosgetsol(Program1,dVa_dth); %Getting solution for dha_dx

temp_surf_axis = [0,4,-2.6,-1.4,0,10];
temp_contour_axis = [-4,4,-3,-1];

figure;
fsurf(SOL_Va,temp_contour_axis,'EdgeColor','none'); % Plot the lyapunov function
%fc = fcontour(SOL_ha, temp_contour_axis,'MeshDensity',200,'Fill','on'); % Plot the lyapunov function
% fc.LevelList = [10 9 8 2 1 0];
% fc.LevelList = [10:0.1:0]
colorbar

%axis(temp_surf_axis+[0,0,0,0,0,300])
xlabel('x');
ylabel('$$\hat{\theta}$$','Interpreter','latex');
zlabel('$$V(x,\hat{\theta})$$','Interpreter','latex')
view(0,90)
saveas(gcf,'barrier_image_positive_values.png')

temp_surf_axis2 = [-4,4,-2.6,-1.4,0,10];


figure;
fsurf(SOL_Va,temp_contour_axis,'EdgeColor','none'); % Plot the lyapunov function
%fc = fcontour(SOL_ha, temp_contour_axis,'MeshDensity',200,'Fill','on'); % Plot the lyapunov function
% fc.LevelList = [10 9 8 2 1 0];
% fc.LevelList = [10:0.1:0]
colorbar
hold on;
plot(PreSet_XTheta_3d,'Color','red','alpha',0.1)

axis(temp_surf_axis)
xlabel('x');
ylabel('$$\hat{\theta}$$','Interpreter','latex');
zlabel('$$V(x,\hat{\theta})$$','Interpreter','latex')
view(0,90)

saveas(gcf,'ha_compared_to_target.png')

figure;
hold on;
[X,Y] = meshgrid(0:0.01:4,-2.6:0.05:-1.4);
ha_eval = double(subs(SOL_Va,struct('x', X, 'theta_h', Y)));
surf(X,Y,ha_eval)
plot(PreSet_XTheta_3d,'Color','red','alpha',0.1)
colorbar

xlabel('x');
ylabel('$$\hat{\theta}$$','Interpreter','latex');
zlabel('$$V(x,\hat{\theta})$$','Interpreter','latex')
view(0,90)

saveas(gcf,'ha_compared_to_target_nofsurf.png')

%% Simulate System

n_sims = 10;
T_sim = 10;
dt = 0.1;
x0 = unifrnd(2,5,[n_sims,1]);
th0 = sampleFromPolytope(scalar_system.Theta,n_sims);
x=[]; u = []; theta = [];

Gamma1 = 1;

for sim_index = 1:n_sims
    x0_i = x0(sim_index);

    x_i = x0_i;
    theta_i = sampleFromPolytope(scalar_system.Theta); 
    u_i = [];

    for k = 1:T_sim-1
        V = double(subs(SOL_Va,struct('x', x_i(k), 'theta_h', theta_i(k) ) ));
        dV_dx = double(subs(SOL_dVa_dx,struct('x', x_i(k), 'theta_h', theta_i(k) ) ));
        dV_dth = double(subs(SOL_dVa_dth,struct('x', x_i(k), 'theta_h', theta_i(k) ) ));
        
        [u_k, opt_out] = aCLF_control( x_i(k) , V, dV_dx, dV_dth , scalar_system , Gamma1 );
        u_i(k) = u_k;
        x_i(:,k+1) = x_i(:,k) + scalar_system.dynamics(x_i(:,k),u_k) * dt;

        G_x = scalar_system.G(x_i(:,k));
        theta_i(:,k+1) = theta_i(:,k)  %+ Gamma1 * dV_dx * (scalar_system.F(x_i(k,:)) + G_x{1} ) * dt ;

    end
    x = [x; x_i];
    u = [u; u_i];
    theta = [theta; theta_i];

end

% Plot the value of the aCLF over time
figure;
hold on;
for sim_index = 1:n_sims % Plot sims
    x_i = x(sim_index,:);
    theta_i = theta(sim_index,:);
    V_i = double(subs(SOL_Va,struct('x', x_i, 'theta_h', theta_i)));
    plot(V_i);
end
xlabel('Simulation Time (k)')
ylabel('$$V_a(x_k)$$','Interpreter','latex')
title('Lyapunov Function Change over time')

saveas(gcf,'simulated_aCLF_Trajectories_const_est.png')


% Plot the trajectories on the surface
figure;
hold on;

% Plot Surfaces
surf(X,Y,ha_eval)
plot(PreSet_XTheta_3d,'Color','red','alpha',0.01)
colorbar

%Plot trajectories on top
for sim_index = 1:n_sims % Plot sims
    x_i = x(sim_index,:);
    theta_i = theta(sim_index,:);
    V_i = double(subs(SOL_Va,struct('x', x_i, 'theta_h', theta_i)));
    plot3(x_i,theta_i,V_i+0.5,':','Color','magenta','LineWidth',1.5);
end

xlabel('x');
ylabel('$$\hat{\theta}$$','Interpreter','latex');
zlabel('$$V(x,\hat{\theta})$$','Interpreter','latex')
view(0,90)

saveas(gcf,'simulated_aCLF_Trajectories_over_surf_const_estimator.png')

%% Create Function For Controller Of This System

function [ u_opt , optim_out ] = aCLF_control( x , V, dV_dx, dV_dth , system_dyn , Gamma )
    %Description:
    %
    %Inputs
    %   V : Current value of the adaptive CLF
    %   dV_dx : Current value of the adaptive CLF's gradient with respect
    %           to the state x
    %   dV_dth: Current value of the adaptive CLF's gradient with respect
    %           to the unknown parameters theta
    %   


    % Constants
    settings_in = sdpsettings('verbose',0);
    dim_x = 1;
    dim_u = 1;

    decay_rate = 0;

    V_Theta = system_dyn.Theta.V;
    n_VTheta = size(V_Theta,1);

    % Algorithm
    u = sdpvar(dim_u,1);
    constraints = [];

    for theta_index = [1:n_VTheta]

        % Collect Vertex
        v_Theta = V_Theta(theta_index,:)';

        % Create constraints

        Lf_V = dV_dx * system_dyn.f(x);
        LF_V = dV_dx * system_dyn.F(x) * v_Theta;
        Lg_V = dV_dx * system_dyn.g(x);
        
        G = system_dyn.G(x);
        LsumG_V = dV_dx*v_Theta(1)*G{1};
        for G_index = [2:length(G)]
            LsumG_V = LsumG_V + dV_dx*v_Theta(G_index)*G{G_index}
        end
        
        % Final constraint
        lyapunov_function_decrease_constraint = [Lf_V + LF_V + (Lg_V + LsumG_V)*u <= - decay_rate * V];
        
        constraints = constraints + lyapunov_function_decrease_constraint;
    end

    % Optimize
    optim_out = optimize(constraints,[ u.^2 ],settings_in);

    u_opt = value(u); % Save the feasible u

end
