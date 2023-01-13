%%  example2_taylor
%   Description:
%       This script visualizes attempts to find an adaptive CBF for the
%       adaptive cruise control system that is discussed in the Taylor and
%       Ames paper "Adaptive Safety with Control Barrier Functions".

clear all; close all; clc;

%% Add Functions to Path
addpath(genpath('../src/'))

%% Constants
acc1 = AdaptiveCruiseControl();

%% Plot Simple Trajectories of the System
x0 = [15;26];
tspan1 = [0,4];
[t1,y1] = ode45( ...
    @(t,x) acc1.dynamics(x,0), ...
    tspan1, x0 ...
    );

[t2,y2] = ode45( ...
    @(t,x) acc1.dynamics(x,-norm(x,2)*0.2), ...
    tspan1, x0 ...
    );

min_x1 = min([y1(:,1);y2(:,1)]);
min_x2 = min([y1(:,2);y2(:,2)]);
max_x1 = max([y1(:,1);y2(:,1)]);
max_x2 = max([y1(:,2);y2(:,2)]);

matching_axis1 = [ tspan1, min_x1, max_x1];

figure;
subplot(1,2,1)
plot(t1,y1(:,1))
xlabel('Time (s)')
ylabel('x')
axis(matching_axis1)
title('Zero Input Trajectory')

subplot(1,2,2)
plot(t2,y2(:,1))
xlabel('Time (s)')
ylabel('x')
axis(matching_axis1)
title('Slight corrective (u = -0.2x) trajectory')

%% Use the CLF-based controller for this work with known theta
%  ===========================================================

% Constants
alpha = 100;
dtheta = [1.0,0.1,0.05];
Theta = Polyhedron('lb',-dtheta,'ub',dtheta) + acc1.theta;

Gamma = eye(3);

deltaV = 0.1;

K_D = 10;
D_star = 30;

K_v = -1;
v_star = acc1.v0;

% Create a dummy program from the tutorial for SOSTools
syms v1 D1 f1 f2 f3 real;

u_min = -9.8*acc1.m;

Program1 = sosprogram([v1;D1]);
monom1 = monomials([v1;D1],[1:4]);
[Program1,l1] = sospolyvar(Program1,monom1);
[Program1,l2] = sospolyvar(Program1,monom1);
ha = alpha.^2 - (D1 - 1.8*v1 - alpha).^2;
Va = (D1 - 1.8*v1).^2;

Verify_aCLF_old([v1;D1],[f1;f2;f3],Va, acc1, Gamma)

return

% %% Use the aCLF-based controller for this work
% %  ===========================================

% syms x a b c d;

% alpha1 = 0.1;

% Program2 = sosprogram([x]);
% %monom1 = monomials([x],[1,2,3,4]);
% [Program2,ha] = sospolyvar(Program2,monom1);

% l1 = 0.01*x^2;
% l2 = 0.01*x^4;

% s1 = x^2;
% [Program2,p2] = sospolyvar(Program2,monom1);

% dV_dx = diff(ha,x);

% % Create constraints
% Program2 = sosineq( Program2 , ha - l1 ); %F2
% tempExpr2_lb = ...
%     - ( dV_dx * ( scalar_system.f(x) + scalar_system.F(x)* Theta_lb ) + ...
%     (dV_dx*( scalar_system.g(x) + scalar_system.G(x) * Theta_lb ) * ( -(1/(Theta_lb+1)) - 1 ) * x ) + l2)
% Program2 = sosineq( Program2 , ...
%     tempExpr2_lb ...
% ); %F3

% tempExpr2_ub = ...
%     - ( dV_dx * ( scalar_system.f(x) + scalar_system.F(x)* Theta_ub ) + ...
%     (dV_dx*( scalar_system.g(x) + scalar_system.G(x) * Theta_ub ) * ( -(1/(Theta_ub+1)) - 1 ) * x ) + l2)
% Program2 = sosineq( Program2 , ...
%     tempExpr2_ub ...
% ); %F3

% % Optimize
% Program2 = sossolve( Program2 );
% assert( (Program1.solinfo.info.pinf ~= 1) || (Program1.solinfo.info.dinf ~= 1) ); % Make sure that neither problem is infeasible?

% % Plot the obtained lyapunov function
% SOL_V2 = sosgetsol(Program2,ha) %Getting solution for V
% figure;
% ezplot(SOL_V2)

% % Plot the lyapunov function over time.
% x0s = [-1.0,-0.6,0.1,1.0];
% thetas = unifrnd(Theta_lb,Theta_ub,1,length(x0s))
% x_h = []; ha_h = [];
% for x0_index = [1:length(x0s)]
%     x0 = x0s(x0_index);
%     scalar_system.theta = thetas(x0_index);

%     x = [x0];
%     ha = [subs(SOL_V,x0)];

%     x_k = x0;
%     ha_k = subs(SOL_V,x_k);
%     dt = 0.1;
%     T  = 10;
%     for k = [0:dt:T-dt]
%         % u_k = ( -(1/(scalar_system.theta+1)) - 1 )*x_k;

%         u_k, optim_out = CLF_control( x , h_at_x, dh_dx_at_x , system_dyn );

%         x_kp1 = x_k + dt * scalar_system.dynamics(x_k,u_k);
%         x = [x,x_kp1];

%         V_kp1 = subs(SOL_V,x_k);
%         ha = [ha,V_kp1];

%         % Update x_k
%         x_k = x_kp1;

%     end

%     x_h = [x_h ; x];
%     ha_h = [ha_h ; ha];
% end

% figure;
% subplot(2,1,1)
% plot([0:dt:T],x_h)
% xlabel('Time')
% ylabel('x')
% title('Unknown Theta, Simple aCLF Control')

% subplot(2,1,2)
% plot([0:dt:T],ha_h)
% xlabel('Time')
% ylabel('V_t')
% title('Unknown Theta, Simple aCLF Control')

% return

% %% Use the aCBF-based controller for this work
% %  ===========================================

% syms v1 D1 a b c d;

% alpha1 = 0.1;

% Program3 = sosprogram([v1 D1]);
% monom1 = monomials([v1 D1],[1,2,3,4]);
% [Program3,V3] = sospolyvar(Program3,monom1);

% l1 = 0.01*x^2;
% l2 = 0.01*x^4;

% s1 = x^2;
% [Program3,p2] = sospolyvar(Program3,monom1);

% dV3_dx = diff(V3,x);

% % Create constraints
% Program3 = sosineq( Program3 , V3 - l1 ); %F2
% tempExpr3_lb = ...
%     - ( dV3_dx * ( scalar_system.f(x) + scalar_system.F(x)* Theta_lb ) + ...
%     (dV3_dx*( scalar_system.g(x) + scalar_system.G(x) * Theta_lb ) * ( -(1/(Theta_lb+1)) - 1 ) * x ) + V3)
% Program3 = sosineq( Program3 , tempExpr3_lb ); %F3

% tempExpr3_ub = ...
%     - ( dV3_dx * ( scalar_system.f(x) + scalar_system.F(x)* Theta_ub ) + ...
%     (dV3_dx*( scalar_system.g(x) + scalar_system.G(x) * Theta_ub ) * ( -(1/(Theta_ub+1)) - 1 ) * x ) + V3)
% Program3 = sosineq( Program3 , ...
%     tempExpr3_ub ...
% ); %F3

% % Optimize
% Program3 = sossolve( Program3 );
% assert( (Program3.solinfo.info.pinf ~= 1) || (Program3.solinfo.info.dinf ~= 1) ); % Make sure that neither problem is infeasible?

% % Plot the obtained lyapunov function
% SOL_V3 = sosgetsol(Program3,V3) %Getting solution for V
% figure;
% ezplot(SOL_V3)

% % Plot the lyapunov function over time.
% x0s = [-1.0,-0.6,0.1,1.0];
% thetas = unifrnd(Theta_lb,Theta_ub,1,length(x0s))
% x_h = []; ha_h = [];
% for x0_index = [1:length(x0s)]
%     x0 = x0s(x0_index);
%     scalar_system.theta = thetas(x0_index);

%     x = [x0];
%     ha = [subs(SOL_V3,x0)];

%     x_k = x0;
%     ha_k = subs(SOL_V3,x_k);
%     dt = 0.1;
%     T  = 10;
%     for k = [0:dt:T-dt]
%         u_k = ( -(1/(scalar_system.theta+1)) - 1 )*x_k;
%         x_kp1 = x_k + dt * scalar_system.dynamics(x_k,u_k);
%         x = [x,x_kp1];

%         V_kp1 = subs(SOL_V3,x_k);
%         ha = [ha,V_kp1];

%         % Update x_k
%         x_k = x_kp1;

%     end

%     x_h = [x_h ; x];
%     ha_h = [ha_h ; ha];
% end

% figure;
% subplot(2,1,1)
% plot([0:dt:T],x_h)
% xlabel('Time')
% ylabel('x')
% title('Unknown Theta, Exponentially Convergent aCLF Control')

% subplot(2,1,2)
% plot([0:dt:T],ha_h)
% xlabel('Time')
% ylabel('V_t')
% title('Unknown Theta, Exponentially Convergent aCLF Control')

function [ u_opt , optim_out ] = CLF_control( x , V_at_x, dV_dx_at_x , system_dyn )

    % Constants
    settings_in = sdpsettings('verbose',0);
    dim_x = 2;
    dim_u = 1;

    D_des = 20; %Desired headway
    decay_rate = 2;

    % Algorithm
    u = sdpvar(dim_u,1);

    Lf_V = dV_dx_at_x * system_dyn.f(x);
    LF_V = dV_dx_at_x * system_dyn.F(x) * system_dyn.theta;
    Lg_V = dV_dx_at_x * system_dyn.g(x) * u;

    lyapunov_function_decrease_constraint = [Lf_V + LF_V + Lg_V <= - decay_rate * V_at_x];

    % Optimize
    constraints = lyapunov_function_decrease_constraint;
    optim_out = optimize(constraints,[ u.^2 ],settings_in);

    u_opt = value(u); % Save the feasible u

end
