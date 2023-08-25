function [V] = BasicACLFSearch(nonlinearSystemWParameters,Gamma)
%FindACLF Uses optimization to attempt to find an aCLF for the origin
%   Detailed explanation goes here

%% Constants


%% Algorithm

% Create YALMIP Optimization Variables

x = sdpvar(1,1);
a = sdpvar(1,1);
l1 = sdpvar(1,1);

V = (x-a)^4;

s1 = x^2;

dV_dx = 4 * x^3;

% Create constraints
% F1 = sos(V); % We know V is an SOS
F2 = sos(V - l1);

outputArg1 = inputArg1;
outputArg2 = Gamma;
end