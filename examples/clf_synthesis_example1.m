% clf_synthesis_example1.m
% Description:
%   Attempting to implement Example 1 in 5.1.4.1 of the Weehong thesis.

clear all; close all; clc;

%% Add Functions to Path
addpath(genpath('../src/'))

%% Constants

% Define System
%   \dot{x} = (Lx)*u

dim_x = 2;
dim_u = 1;

L = [ 3, 5; -20, 10 ];

%% Create SOS Program Using YALMIP

% x = sdpvar(2,1);
% P = sdpvar(2,2,'symmetric');
% 
% % Create first constraint
% [V,c_V,v_V] = polynomial(x,2);
% l1 = (1e-4)*x'*x;
% 
% constr15_expr = sos(V-l1); %Named 15 after (5.15) in Weehong's thesis
% 
% % Create second constraint
% [s1,c_s1,v_s1] = polynomial(x,4);
% [p2,c_p2,v_p2] = polynomial(x,2);
% 
% dV_dx = jacobian(V,x);
% 
% f_x = zeros(dim_x,1);
% 
% l2 = (0.01)*x'*x;
% constr16_expr = sos(- ( s1 * ( dV_dx * f_x ) + p2 * (dV_dx*L*x) + l2 ));
% 
% % Optimize
% options = sdpsettings('sos.newton',0,'sos.congruence',0,'debug',1);
% [sol,v,Q] = solvesos([constr15_expr,constr16_expr,sos(s1)],[],options,[c_s1;c_p2;c_V])

% x = sdpvar(1,1);y = sdpvar(1,1);
% p = (1+x)^4 + (1-y)^2;
% F = sos(p);
% [sol,v,Q] = solvesos(F);
% 
% 
% sdpvar x y t
% p1 = t*(1+x*y)^2-x*y+(1-y)^2;
% p2 = (1-x*y)^2+x*y+t*(1+y)^2;
% F = [sos(p1), sos(p2)];
% solvesos(F,t);
% value(t)

%% Create SOS Program Using SOSTOOLS

% Create a dummy program from the tutorial for SOSTools
syms x1 x2 real;

Program1 = sosprogram([x1;x2]);
monom_V = monomials([x1;x2],[1:2]);
[Program1,V] = sospolyvar(Program1,monom_V);

l1 = 0.1*(x1.^2 + x2.^2);

monom_l2 = monomials([x1;x2],[0:4]);
[Program1,l2] = sossosvar(Program1,monom_l2);
%l2 = 0.01*(x1.^4 + x2.^4);

Program1 = sosineq( Program1 , V - l1 ); 

% Second Constraint
monom_s1 = monomials([x1;x2],[0:4]);
[Program1,s1] = sospolyvar(Program1,monom_s1);

monom_p2 = monomials([x1;x2],[0:2]);
[Program1,p2] = sospolyvar(Program1,monom_p2);

dV_dx(1,1) = diff(V,x1);
dV_dx(1,2) = diff(V,x2);

f_x = zeros(dim_x,1);

p2 * (dV_dx*L*[x1;x2])
%constr16_expr = - (  s1 * (dV_dx * f_x ) + (x1.^2 + x2.^2) * (dV_dx*L*[x1;x2]) + l2 );
constr16_expr = - ( (x1.^2 + x2.^2)*(dV_dx*L*[x1;x2]) - l2 );
Program1 = sosineq( Program1 , constr16_expr );

% Optimize
solver_opt.solver = 'sedumi';
Program1 = sossolve( Program1 , solver_opt );
assert( (Program1.solinfo.info.pinf ~= 1) || (Program1.solinfo.info.dinf ~= 1) ); % Make sure that neither problem is infeasible?
SOL_V = sosgetsol(Program1,V);

figure;
ezcontour(SOL_V)

figure;
ezsurf(SOL_V)