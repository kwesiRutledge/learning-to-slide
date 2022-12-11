%% max_ellipsoid_test.m
%   Description:
%       This function should create a suite of tests for the ellipsoid
%       algorithm.

clear all; close all; clc;

% Constants

A = [ eye(2) ; -eye(2) ];
b = [1.4,1.3,1.2,1.1]';

X_target = Polyhedron(A,b)
Theta = Polyhedron('lb',0.5,'ub',0.8)

% Algorithm

P = sdpvar(3,3);
c = sdpvar(3,1,'full');

constraints = [ P >= 0 ];

for row_index = [1:size(A,1)]
    %Extract row from A and b
    A_i = A(row_index,:);
    b_i = b(row_index);

    A_i_tilde = [ A_i , zeros(1,Theta.Dim) ];
    constraints = constraints + [ A_i_tilde * P * P' * A_i_tilde' <= b_i - A_i_tilde*c ];
end

for row_index = [1:size(Theta.A,1)]
    %Extract row from A and b
    A_i = Theta.A(row_index,:);
    b_i = Theta.b(row_index);

    A_i_tilde = [ zeros(1,size(X_target.A,2)) , A_i ];
    constraints = constraints + [ A_i_tilde * P * P' * A_i_tilde' <= (b_i - A_i_tilde*c)'*(b_i - A_i_tilde*c) ];
end

% for vertex_index = 1:size(Theta.V,1)
%     % Add in the constraints on containment of the 
% end

output__struct = optimize( ...
    constraints, ...
    -logdet(P) ...
    )

% Test 2

Theta = Polyhedron('lb',-0.5,'ub',0.8)

P2 = sdpvar(3,3);
% c = sdpvar(3,1,'full');

constraints = [ P >= 0 ];

for row_index = [1:size(A,1)]
    %Extract row from A and b
    A_i = A(row_index,:);
    b_i = b(row_index);

    A_i_tilde = [ A_i , zeros(1,Theta.Dim) ];
    constraints = constraints + [ A_i_tilde * P2 * A_i_tilde' <= b_i.^2 ];
end

for row_index = [1:size(Theta.A,1)]
    %Extract row from A and b
    A_i = Theta.A(row_index,:);
    b_i = Theta.b(row_index);

    A_i_tilde = [ zeros(1,size(X_target.A,2)) , A_i ];
    constraints = constraints + [ norm(A_i_tilde*P2,2) <= b_i.^2 ];
end

% for vertex_index = 1:size(Theta.V,1)
%     % Add in the constraints on containment of the 
% end

output__struct = optimize( ...
    constraints, ...
    -logdet(P) ...
    )