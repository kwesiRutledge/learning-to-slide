%% max_ellipsoid_test.m
%   Description:
%       This function should create a suite of tests for the ellipsoid
%       algorithm.

clear all; close all; clc;

%% Add Library Files
addpath(genpath('../../src/'))


%% Constants

A = [ eye(2) ; -eye(2) ];
b = [1.4,1.3,1.2,1.1]';

X_target = Polyhedron(A,b);
Theta = Polyhedron('lb',0.5,'ub',0.8);
system1 = SimpleSystem1(Theta);
assert( Theta.contains( system1.theta ) )

% New Settings
ops = sdpsettings('debug','1');

% Tests
% sdpvar x y t
% p1 = t*(1+x*y)^2-x*y+(1-y)^2;
% p2 = (1-x*y)^2+x*y+t*(1+y)^2;
% F = [sos(p1), sos(p2)];
% solvesos(F,t);
% value(t)

%% Algorithm

x = sdpvar(1,1);
%a = sdpvar(1,1); b = sdpvar(1,1); c = sdpvar(1,1); d = sdpvar(1,1);

l1 = 0.01*x^2;
l2 = 0.01*x^2;

[V, c1] = polynomial(x,4)

s1 = x^2;
% p2 = a + b*x + c*x^2 + d*x^3;
[p2,c2,v] = polynomial([x],4);

sdisplay(p2)

derivative_mat = zeros(length(v));
for dm_index = 1:size(derivative_mat,1)-1
    derivative_mat(dm_index,dm_index+1) = dm_index;
end

c3 = derivative_mat * c2;
dV_dx = c3'*v;

% Create constraints
% F1 = sos(V); % We know V is an SOS
F2 = sos(V - l1);
F3 = sos( ...
    - ( s1 * dV_dx * ( system1.f(x) + system1.F(x)* system1.theta ) + p2 * (dV_dx*( system1.g(x) + system1.G(x) * system1.theta ) ) + l2) ...
    );
F = [F2,F3];

% Call Solver
[sol,v,Q] = solvesos(F,[],ops,[c1,c2])