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
        U;
        X0;

        dim_x;
        dim_u;
        dim_theta;
    end
    
    methods
        function obj = SimpleSystem1(varargin)
            %SIMPLESYSTEM1 Construct an instance of this class
            %   Provide the unknown parameter set
            %   Polyhedron() object
            %Usage:
            %   system1 = SimpleSystem('Theta',Polyhedron('lb',0.5,'ub',0.8));
            %   system1 = SimpleSystem('Theta',Polyhedron('lb',0.5,'ub',0.8),'theta',0.7)
            %   system1 = SimpleSystem('Theta',Theta,'theta',theta)
            
            [ Theta , U , X0, ss_settings ] = input_processing_SimpleSystem1(varargin{:});

            
            obj.Theta = Theta;
            obj.theta = ss_settings.theta;

            obj.U = U;
            obj.X0 = X0;
            % Dimensions of State and input
            obj.dim_x = 1;
            obj.dim_u = 1;
            obj.dim_theta = 1;

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

function [Theta, U , X0, ss_settings] = input_processing_SimpleSystem1(varargin)
    %Description:
    %   Process the inputs given to ExternalBehaviorSet() constructor.

    %% Set Defaults

    ss_settings = struct( ...
        'Theta', Polyhedron('lb',0.5,'ub',0.8), ...
        'theta', sampleFromPolytope(Polyhedron('lb',0.5,'ub',0.8)), ...
        'U', Polyhedron('lb',-10,'ub',10), ...
        'X0', Polyhedron('lb',1.0,'ub',2.0) ...
    );

    %% Algorithm
    if nargin > 0
        v_index = 1;
        switch varargin{v_index}
            case 'Theta'
                ss_settings.Theta = varargin{v_index+1};
                v_index = v_index + 2;
            case 'theta'
                ss_settings.theta = varargin{v_index+1};
                if ~ss_settings.Theta.contains(ss_settings.theta)
                    error(['The theta input that was given is not inside the set ''Theta''. Please propose a new one'])
                end
                v_index = v_index + 2;
            otherwise
                error(['Unexpected input to SimpleSystem1: ' varargin{v_index} ])
        end

    end

    %% Create Outputs
    Theta = ss_settings.Theta;
    U = ss_settings.U;
    X0 = ss_settings.X0;
end

