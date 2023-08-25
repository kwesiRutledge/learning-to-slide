%%  example1_scalar
%   Description:
%       This script visualizes some useful scripts regarding the simple
%       system example developed for our MATLAB implementation.

clear all; close all; clc;

%% Add Functions to Path
addpath(genpath('../../src/'))

%% Constants
Theta_lb = 0.5; Theta_ub = 0.85;
Theta1 = Polyhedron('lb',Theta_lb,'ub',Theta_ub);
scalar_system = SimpleSystem1('Theta', Theta1);

%% Plot Simple Trajectories of the System
tspan1 = [0,4];
[t1,y1] = ode45( ...
    @(t,x) scalar_system.dynamics(x,0), ...
    tspan1, 0.1 ...
    );

[t2,y2] = ode45( ...
    @(t,x) scalar_system.dynamics(x,-x*0.2), ...
    tspan1, 0.1 ...
    );

matching_axis = [ tspan1, min([y1;y2]), max([y1;y2])];

figure;
subplot(1,2,1)
plot(t1,y1)
xlabel('Time (s)')
ylabel('x')
axis(matching_axis)
title('Zero Input Trajectory')

subplot(1,2,2)
plot(t2,y2)
xlabel('Time (s)')
ylabel('x')
axis(matching_axis)
title('Slight corrective (u = -0.2x) trajectory')

%% Use the CLF-based controller for this work with known theta
%  ===========================================================

% Create a dummy program from the tutorial for SOSTools
syms x a b c d;

Program1 = sosprogram([x]);
monom1 = monomials([x],[1,2,3,4]);
[Program1,V] = sospolyvar(Program1,monom1);

l1 = 0.01*x^2;
l2 = 0.01*x^4;

s1 = x^2;
[Program1,p2] = sospolyvar(Program1,monom1);

dV_dx = diff(V,x);

% Create constraints
Program1 = sosineq( Program1 , V - l1 ); %F2
tempExpr = ...
    - ( dV_dx * ( scalar_system.f(x) + scalar_system.F(x)* scalar_system.theta ) + ...
    (dV_dx*( scalar_system.g(x) + scalar_system.G(x) * scalar_system.theta ) * ( -(1/(scalar_system.theta+1)) - 1 ) * x ) + l2)
Program1 = sosineq( Program1 , ...
    tempExpr ...
); %F3

% Optimize
Program1 = sossolve( Program1 );
assert( (Program1.solinfo.info.pinf ~= 1) || (Program1.solinfo.info.dinf ~= 1) ); % Make sure that neither problem is infeasible?

% Plot the obtained lyapunov function
SOL_V = sosgetsol(Program1,V); %Getting solution for V
figure;
ezplot(SOL_V) % Plot the lyapunov function
ylabel('$$V(x)$$','Interpreter','latex')
saveas(gcf,'test1-V.png')

% Plot the lyapunov function over time.
x0s = [-1.0,-0.6,0.1,1.0];
x_h = []; V_h = [];
for x0_index = [1:length(x0s)]
    x0 = x0s(x0_index);
    x = [x0];
    V = [subs(SOL_V,x0)];

    x_k = x0;
    V_k = subs(SOL_V,x_k);
    dt = 0.1;
    T  = 10;
    for k = [0:dt:T-dt]
        u_k = ( -(1/(scalar_system.theta+1)) - 1 )*x_k;
        x_kp1 = x_k + dt * scalar_system.dynamics(x_k,u_k);
        x = [x,x_kp1];

        V_k = subs(SOL_V,x_k);
        V = [V,V_k];

        % Update x_k
        x_k = x_kp1;

    end

    x_h = [x_h ; x];
    V_h = [V_h ; V];
end

figure;
subplot(2,1,1)
plot([0:dt:T],x_h)
xlabel('Time')
ylabel('x')

subplot(2,1,2)
plot([0:dt:T],V_h)
xlabel('Time')
ylabel('V_t')
saveas(gcf,'test1-simulation.png')

%% Use the aCLF-based controller for this work
%  ===========================================

syms x a b c d;

alpha1 = 0.1;

Program2 = sosprogram([x]);
%monom1 = monomials([x],[1,2,3,4]);
[Program2,V] = sospolyvar(Program2,monom1);

l1 = 0.01*x^2;
l2 = 0.01*x^4;

s1 = x^2;
[Program2,p2] = sospolyvar(Program2,monom1);

dV_dx = diff(V,x);

% Create constraints
Program2 = sosineq( Program2 , V - l1 ); %F2
tempExpr2_lb = ...
    - ( dV_dx * ( scalar_system.f(x) + scalar_system.F(x)* Theta_lb ) + ...
    (dV_dx*( scalar_system.g(x) + scalar_system.G(x) * Theta_lb ) * ( -(1/(Theta_lb+1)) - 1 ) * x ) + l2);
Program2 = sosineq( Program2 , ...
    tempExpr2_lb ...
); %F3

tempExpr2_ub = ...
    - ( dV_dx * ( scalar_system.f(x) + scalar_system.F(x)* Theta_ub ) + ...
    (dV_dx*( scalar_system.g(x) + scalar_system.G(x) * Theta_ub ) * ( -(1/(Theta_ub+1)) - 1 ) * x ) + l2);
Program2 = sosineq( Program2 , ...
    tempExpr2_ub ...
); %F3

% Optimize
Program2 = sossolve( Program2 );
assert( (Program1.solinfo.info.pinf ~= 1) || (Program1.solinfo.info.dinf ~= 1) ); % Make sure that neither problem is infeasible?

% Plot the obtained lyapunov function
SOL_V2 = sosgetsol(Program2,V) %Getting solution for V
figure;
ezplot(SOL_V2)
ylabel('V(x)')
saveas(gcf,'test2-V.png')

% Plot the lyapunov function over time.
x0s = [-1.0,-0.6,0.1,1.0];
thetas = unifrnd(Theta_lb,Theta_ub,1,length(x0s))
x_h = []; V_h = [];
for x0_index = [1:length(x0s)]
    x0 = x0s(x0_index);
    scalar_system.theta = thetas(x0_index);

    x = [x0];
    V = [subs(SOL_V,x0)];

    x_k = x0;
    V_k = subs(SOL_V,x_k);
    dt = 0.1;
    T  = 10;
    for k = [0:dt:T-dt]
        u_k = ( -(1/(scalar_system.theta+1)) - 1 )*x_k;
        x_kp1 = x_k + dt * scalar_system.dynamics(x_k,u_k);
        x = [x,x_kp1];

        V_kp1 = subs(SOL_V,x_k);
        V = [V,V_kp1];

        % Update x_k
        x_k = x_kp1;

    end

    x_h = [x_h ; x];
    V_h = [V_h ; V];
end

figure;
subplot(2,1,1)
plot([0:dt:T],x_h)
xlabel('Time')
ylabel('x')
title('Unknown Theta, Simple aCLF Control')

subplot(2,1,2)
plot([0:dt:T],V_h)
xlabel('Time')
ylabel('V_t')
title('Unknown Theta, Simple aCLF Control')
saveas(gcf,'test2-simulation.png')

%% Use the aCLF-based controller for this work
%  ===========================================

syms x theta a b c d;

alpha1 = 0.1;
Gamma3 = eye(1);

Program3 = sosprogram([x;theta]);
% monom3 = monomials([x;theta],[1,2,3,4]);
% [Program3,V3] = sospolyvar(Program3,monom3);
V3 = 0.6*x^2 + 0.01*x^3 + 0.3296 * x^4 + 0.1 * theta^4;

monom_l_i = monomials([x; theta],[1,2,3,4]);
[Program3,l1] = sospolyvar(Program3,monom_l_i);
%l1 = 0.01*(x^2+theta^2);
[Program3,l2] = sospolyvar(Program3,monom_l_i);
%l2 = 0.01*x^4+(1e-4)*theta^4;

[Program3,alpha2] = sospolyvar(Program3,monom_l_i);

s1 = x^2;
[Program3,p2] = sospolyvar(Program3,monom1);

dV3_dx = diff(V3,x);
dV3_dth = diff(V3,theta);

% Create constraints
Program3 = sosineq( Program3 , V3 - l1 ); %F2
tempExpr3_lb = ...
    - ( dV3_dx * ( scalar_system.f(x) + scalar_system.F(x) * (Theta_lb + Gamma3 * dV3_dth ) ) + ...
    (dV3_dx*( scalar_system.g(x) + scalar_system.G(x) * (Theta_lb + Gamma3 * dV3_dth ) ) * ( -(1/(Theta_lb+1)) - 1 ) * x ) + alpha2*V3);
Program3 = sosineq( Program3 , tempExpr3_lb ); %F3

tempExpr3_ub = ...
    - ( dV3_dx * ( scalar_system.f(x) + scalar_system.F(x) * (Theta_ub + Gamma3 * dV3_dth ) ) + ...
    (dV3_dx*( scalar_system.g(x) + scalar_system.G(x) * (Theta_ub + Gamma3 * dV3_dth ) ) * ( -(1/(Theta_ub+1)) - 1 ) * x ) + alpha2 * V3);
Program3 = sosineq( Program3 , ...
    tempExpr3_ub ...
); %F3

% Optimize
Program3 = sossolve( Program3 );
assert( (Program3.solinfo.info.pinf ~= 1) && (Program3.solinfo.info.dinf ~= 1) ); % Make sure that neither problem is infeasible?

% Plot the obtained lyapunov function
SOL_V3 = sosgetsol(Program3,V3) %Getting solution for V
figure;
ezsurf(SOL_V3)
ylabel('$$V(x)$$','Interpreter','latex')
saveas(gcf,'test3-V.png')

% Plot the lyapunov function over time.
x0s = [-1.0,-0.6,0.1,1.0];
thetas = unifrnd(Theta_lb,Theta_ub,1,length(x0s))
x_h = []; V_h = [];
for x0_index = [1:length(x0s)]
    x0 = x0s(x0_index);
    scalar_system.theta = thetas(x0_index);

    x = [x0];
    V = [subs(SOL_V3,x0)];

    x_k = x0;
    V_k = subs(SOL_V3,x_k);
    dt = 0.1;
    T  = 10;
    for k = [0:dt:T-dt]
        u_k = ( -(1/(scalar_system.theta+1)) - 1 )*x_k;
        x_kp1 = x_k + dt * scalar_system.dynamics(x_k,u_k);
        x = [x,x_kp1];

        V_kp1 = subs(SOL_V3,x_k);
        V = [V,V_kp1];

        % Update x_k
        x_k = x_kp1;

    end

    x_h = [x_h ; x];
    V_h = [V_h ; V];
end

figure;
subplot(2,1,1)
plot([0:dt:T],x_h)
xlabel('Time')
ylabel('x')
title('Unknown Theta, Exponentially Convergent aCLF Control')

subplot(2,1,2)
plot([0:dt:T],V_h)
xlabel('Time')
ylabel('V_t')
title('Unknown Theta, Exponentially Convergent aCLF Control')
saveas(gcf,'test3-simulation.png')

function [ u_opt , optim_out ] = aCLF_control( x , V_at_x, dV_dx_at_x , system_dyn )

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
