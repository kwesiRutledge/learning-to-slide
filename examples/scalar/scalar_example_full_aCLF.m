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
scalar_system = SimpleSystem1('Theta',Theta1);
scalar_system.U = scalar_system.U * 2;

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

%% Use the aCLF-based controller for this work
%  ===========================================

syms x theta a b c d;

alpha1 = 0.1;
Gamma3 = eye(1);

Program3 = sosprogram([x;theta]);

monom1 = monomials([x],[1,2,3,4]);

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
zlabel('$$V(x)$$','Interpreter','latex')
saveas(gcf,'images/scalar-capa2-full-aclf-V.png')

%% Simulate the System and Plot The Controller's performance

T_sim = 30;
N_sims = 2;
decay_rate = 1;
X0 = Polyhedron('lb',1.0,'ub',2.0);
dt = 0.1;

% Plot the lyapunov function over time.
[x, th0, V_history, th_history, th_star, u_history] = simulate_capa2( ...
    scalar_system, ...
    SOL_V3, ...
    T_sim, X0 , dt, ...
    Gamma3, ...
    N_sims, decay_rate, ...
    {'x'}, {'theta'});

figure;
subplot(2,1,1)
hold on;
for sim_index = 1:N_sims
    plot([0:dt:T_sim*dt],x{sim_index})
end
xlabel('Time')
ylabel('x')
title('Unknown Theta, Exponentially Convergent aCLF Control')

subplot(2,1,2)
hold on;
for sim_index = 1:N_sims
    plot([0:dt:T_sim*dt],th_history{sim_index})
end
saveas(gcf,'images/scalar-capa2-x-theta-simulation.png')

figure;
hold on;
for sim_index = 1:N_sims
    plot([0:dt:T_sim*dt],V_history{sim_index})
end
xlabel('Time')
ylabel('V_t')
title('Unknown Theta, Exponentially Convergent aCLF Control')
saveas(gcf,'images/scalar-capa2-Va-simulation.png')


%% Save data
date_string = datestr(now,'ddmmmyyyy-HHMM');
save(['data/scalar-capa2-full-aCLF-' date_string '.mat' ])
