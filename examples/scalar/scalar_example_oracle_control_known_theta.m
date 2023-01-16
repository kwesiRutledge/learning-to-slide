%%  scalar_example_oracle_control_known_theta
%   Description:
%       This script attempts to find a Lyapunov function for the scalar
%       CAPA2 system ASSUMING that:
%       - parameter value is known
%       - a simple controller is given for the system assuming theta is
%       known.

clear all; close all; clc;

%% Add Functions to Path
addpath(genpath('../../src/'))

%% Constants
Theta_lb = 0.5; Theta_ub = 0.85;
Theta1 = Polyhedron('lb',Theta_lb,'ub',Theta_ub);
scalar_system = SimpleSystem1(Theta1);

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
G = scalar_system.G(x);
G1 = G{1};
tempExpr = ...
    - ( dV_dx * ( scalar_system.f(x) + scalar_system.F(x)* scalar_system.theta ) + ...
    (dV_dx*( scalar_system.g(x) + G1 * scalar_system.theta ) * ( -(1/(scalar_system.theta+1)) - 1 ) * x ) + l2)
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
saveas(gcf,'images/oracle-control-knownTheta-V.png')

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
saveas(gcf,'images/oracle-control-knownTheta-simulation.png')

date_string = datestr(now,'ddmmmyyyy-HHMM');
save(['data/scalar-capa2-knownTheta-oracleControl-' date_string '.mat' ])
