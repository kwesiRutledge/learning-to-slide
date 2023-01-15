function [samples_out] = sampleFromPolytope(varargin)
%SAMPLEFROMPOLYTOPE Summary of this function goes here
%   Samples n possible elements of the set represented by P using
%   the uniform distribution.
%Usage:
%   singleSample = sampleFromPolytope(P1);
%   NSamples = sampleFromPolytope(P1,n);
%
%Outputs
%   samples_out: An n x n_samples matrix which defines n_samples samples
%   from the polytope P1. n is the dimension of the Polytope P1.

%% Input Processing
P = varargin{1};
if ~isa(P,'Polyhedron')
    error('The first input to sampleFromPolytope is not a Polyhedron object!')
end

if ~P.isBounded()
    error('The Polyhedron should be bounded to get a proper result.')
end

n_samples = 1; %Default value

if nargin == 2
    n_samples = varargin{2};
    if ~isnumeric(n_samples)
        error('The second argument to sampleFromPolytope must be numeric!')
    end
end

if nargin > 2
    error('There should be 2 inputs to sampleFromPolyhedron at max!')
end

%% Constants
V = P.V;
n_verts = size(V,1);
mu0 = 1;

%% Algorithm
if n_verts == 1
    %If there is only one vertex given, then return it.
    %It is the only thing in this set.
    sample_out = polytope_in.V';
    return
end

convex_comb = exprnd(mu0,n_verts,n_samples);
convex_comb = convex_comb./repmat(sum(convex_comb),n_verts,1);

samples_out = P.V'*convex_comb;


end

