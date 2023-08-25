%clf_synthesis_example2.m
% Description:
%   Attempting to implement Example 1 in 5.1.4.1 of the Weehong thesis using YALMIP.

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

x = sdpvar(2,1);
P = sdpvar(2,2,'symmetric');

% Create first constraint
% [V,c_V,v_V] = polynomial(x,2);
V = 3.01*x(1).^2 - 0.143*x(1)*x(2) + 1 * x(2).^2;
l1 = (1e-4)*x'*x;

constr15_expr = sos(V-l1); %Named 15 after (5.15) in Weehong's thesis

% Create second constraint
[s1,c_s1,v_s1] = polynomial(x,4);
[p2,c_p2,v_p2] = polynomial(x,2);

dV_dx = jacobian(V,x);

f_x = zeros(dim_x,1);

l2 = (0.01)*x'*x;
constr16_expr = sos(- ( s1 * ( dV_dx * f_x ) + p2 * (dV_dx*L*x) + l2 ));

% Optimize
options = sdpsettings('sos.newton',0,'sos.congruence',0);
options = sdpsettings(options,'debug',1);
[sol,v,Q] = solvesos([constr15_expr,constr16_expr,sos(s1)],[],options,[c_s1;c_p2])

% x = sdpvar(1,1);y = sdpvar(1,1);
% p = (1+x)^4 + (1-y)^2;
% F = sos(p);
% [sol,v,Q] = solvesos(F);


% sdpvar x y t
% p1 = t*(1+x*y)^2-x*y+(1-y)^2;
% p2 = (1-x*y)^2+x*y+t*(1+y)^2;
% F = [sos(p1), sos(p2)];
% solvesos(F,t);
% value(t)
