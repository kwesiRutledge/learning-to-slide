classdef SimpleSystem1
    %SIMPLESYSTEM1 A scalar form of a nonlinear system with unknown
    %parameters that we will study for theory purposes.
    %   This creates a model of the system:
    %       \dot{x} = (1 + a) x + (1 + a)u
    %   and provides a straight-forward method for simulating it.
    %   In this case, the parameter theta is a scalar and the two functions
    %   F(x) and G(x) are identical.
    %   We decompose this model into the form
    %       \dot{x} = f(x) + F(x) \theta + ( g(x) + \sum_i \theta_i G_i(x)) u
    %   from our paper
    
    properties
        Theta;
        theta;
        x;
    end
    
    methods
        function obj = SimpleSystem1(varargin)
            %SIMPLESYSTEM1 Construct an instance of this class
            %   Provide the unknown parameter set
            %   Polyhedron() object
            %Usage:
            %   system1 = SimpleSystem(Polyhedron('lb',0.5,'ub',0.8));
            %   system1 = SimpleSystem(Polyhedron('lb',0.5,'ub',0.8),'theta',0.7)
            
            % Input Checking
            if nargin < 1
                error(['Require at least one input to SimpleSystem! Received ' num2str(nargin)])
            end

            % Extract necessary arguments
            obj.Theta = varargin{1};

            % Assign defaults
            convex_comb_vector = exprnd(1,length(obj.Theta.V),1);
            convex_comb_vector = convex_comb_vector/sum(convex_comb_vector);
            obj.theta = obj.Theta.V' * convex_comb_vector;


            % If there is more than one input, then accept the other
            % values.
            arg_index = 2;
            if nargin > 1
                switch varargin{arg_index}
                    case 'theta'
                        obj.theta = varargin{arg_index+1};
                        arg_index = arg_index + 2;
                    otherwise
                        error(['Unexpected input: ' varargin{arg_index} ])
                end
            end

        end
        
        function f_x = f(obj,x)
            %f() The function represents the component f(x) in our problem
            %structure.
            %   Detailed explanation goes here
            f_x = x;
        end

        function F_x = F(obj,x)
            %F() The function represents the component F(x) in our problem
            %structure.
            %   Detailed explanation goes here
            F_x = x;
        end

        function g_x = g(obj,x)
            %g() The function represents the component g(x) in our problem
            %structure.
            %   Detailed explanation goes here
            g_x = 1;
        end

        function G_x = G(obj,x)
            %G() The function represents the set of all G_i(x) in our problem
            %structure.
            %   For the scalar system, there is only one of these.

            G_x = {1};
        end

        function dx_dt = dynamics(obj,x,u)
            %dynamics The function that computes the overall derivative of
            %the state vector (really a scalar).

            % Constants
            theta = obj.theta;

            % Collect G
            G = obj.G(x);
            sum_Gi = theta(1) * G{1};
            for theta_dim = [2:length(theta)]
                sum_Gi = sum_Gi + theta(theta_dim) * G{theta_dim};
            end

            % Compute derivative
            dx_dt = obj.f(x) + obj.F(x) * theta + ...
                ( obj.g(x) + sum_Gi )* u;
        end
    end
end

