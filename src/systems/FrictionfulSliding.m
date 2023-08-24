classdef FrictionfulSliding
    %FRICTIONFULSLIDING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Spaces
        Theta;  % Parameter space
        X0;     % Initial State space
        U;      % Input space
        X_T;     % Target Set
        % Unkown Parameter
        theta;
        % Known parameters
        m;
        % Theoretical quantities
        dim_x;
        dim_u;
        dim_theta;
    end
    
    methods
        function obj = FrictionfulSliding(varargin)
            %FRICTIONFULSLIDING Construct an instance of this class
            %Usage:
            %   system1 = SimpleSystem('Theta',Polyhedron('lb',0.5,'ub',0.8));
            %   system1 = SimpleSystem('Theta',Polyhedron('lb',0.5,'ub',0.8),'theta',0.7)
            %   system1 = SimpleSystem('Theta',Theta,'theta',theta)
            
            [ Theta , U , X0, X_T, fs_settings ] = input_processing_FrictionfulSliding(varargin{:});

            
            obj.Theta = Theta;
            obj.theta = fs_settings.theta;

            obj.U = U;
            obj.X0 = X0;
            obj.X_T = X_T;

            obj.m = fs_settings.m;

            % Dimensions of State and input
            obj.dim_x = 1;
            obj.dim_u = 1;
            obj.dim_theta = 1;
        end
        
        function f_x = f(obj,x)
            %f() The function represents the component f(x) in our problem
            %structure.
            %   Inputs
            %       x - 2 dimensional state
            
            % Constants
            p_x = x(1);
            v_x = x(2);

            % Algorithm
            f_x = [v_x; 0];
        end

        function F_x = F(obj,x)
            %F() The function represents the component F(x) in our problem
            %structure.
            %   Inputs
            %       x - 2 dimensional state
            
            % Constants
            p_x = x(1);
            v_x = x(2);
            g = 10; %m/s^2

            % Algorithm
            if v_x > 0.01
                F_x = [0;-g * sign(v_x)];
            else
                F_x = zeros(2,1);
            end
        end

        function g_x = g(obj,x)
            %g() The function represents the component g(x) in our problem
            %structure.
            %   Inputs
            %       x - 2 dimensional state
            
            % Constants
            p_x = x(1);
            v_x = x(2);
            m = obj.m;

            % Algorithm
            g_x = [0,0;(1/m),0];
        end

        function G_x = G(obj,x)
            %G() The function represents the set of all G_i(x) in our problem
            %structure.
            %   For the scalar system, there is only one of these.
            %   Inputs
            %       x - 2 dimensional state
            
            % Constants
            p_x = x(1);
            v_x = x(2);

            m = obj.m;

            % Algorithm
            G_x = {[0,0;0,-(1/m)]};
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

        function [plot_handles] = show(obj,x,object_color)
            %Description:
            %   This function shows the sliding object system as we
            %   envisioned it in the research.

            % Constants
            sliding_block_height = 0.5;

            X0_lb = min(obj.X0.V(:,1));
            XT_ub = max(obj.X_T.V(:,1));
            min_table_width = XT_ub - X0_lb;
            desired_axis = [ ...
                X0_lb-0.2*min_table_width ...
                XT_ub+0.2*min_table_width ...
                -0.5, sliding_block_height*2.0];

            X0_color = 'blue';
            X0_linewidth = 2.0;
            X0_line_height = 0.5*sliding_block_height;

            table_linewidth = 4.0;

            X_T_color = 'red';
            X_T_linewidth = X0_linewidth;
            X_T_line_height = X0_line_height;

            if ~exist('object_color')
                object_color = 'black';
            end
            object_lw = 3.0;

            % Draw Initial Condition Set
            plot_handles(1) = plot(X0_lb*ones(2,1),[0,X0_line_height], ...
                'LineStyle',':','Color', X0_color, ...
                'LineWidth',X0_linewidth);
            hold on;
            plot_handles(2) = plot(max(obj.X0.V(:,1))*ones(2,1),[0,X0_line_height], ...
                'LineStyle',':','Color', X0_color, ...
                'LineWidth',X0_linewidth);

            % Draw Target Condition Set
            plot_handles(3) = plot(XT_ub*ones(2,1),[0,X_T_line_height], ...
                'LineStyle',':','Color', X_T_color, ...
                'LineWidth',X_T_linewidth);
            plot_handles(4) = plot(min(obj.X_T.V(:,1))*ones(2,1),[0,X_T_line_height], ...
                'LineStyle',':','Color', X_T_color, ...
                'LineWidth',X_T_linewidth);

            % Draw Table
            plot_handles(5) = plot([desired_axis(1),desired_axis(2)],zeros(1,2), ...
                'LineWidth',table_linewidth);

            % Draw Object
            corners = [ 
                -0.5*sliding_block_height,0;
                -0.5*sliding_block_height,sliding_block_height;
                0.5*sliding_block_height,sliding_block_height;
                0.5*sliding_block_height,0
                ];
            corner_offset = zeros(4,2);
            corner_offset(:,1) = x(1);

            shifted_corners = corners + corner_offset;

            for corner_index = [1:3]
                plot_handles(5+corner_index) = plot( ...
                    shifted_corners([corner_index:corner_index+1],1), ...
                    shifted_corners([corner_index:corner_index+1],2), ...
                    'Color',object_color, ...
                    'LineStyle','-', ...
                    'LineWidth',object_lw);
            end
            plot_handles(9) = plot( ...
                    shifted_corners([4,1],1), shifted_corners([4,1],2), ...
                    'Color',object_color, ...
                    'LineStyle','-', ...
                    'LineWidth',object_lw);

            % Force axis
            axis(desired_axis)

        end
    end
end

function [Theta, U , X0, X_T, fs_settings] = input_processing_FrictionfulSliding(varargin)
    %Description:
    %   Process the inputs given to FrictionfulSliding() constructor.

    %% Set Defaults

    fs_settings = struct( ...
        'Theta', Polyhedron('lb',0.5,'ub',0.8), ...
        'theta', sampleFromPolytope(Polyhedron('lb',0.5,'ub',0.8)), ...
        'U', Polyhedron('lb',-10*ones(1,2),'ub',10*ones(1,2)), ...
        'X0', Polyhedron('lb',[-0.2,1.2],'ub',[0.2,1.8]), ...
        'm', 1.0, ... %kg
        'X_T', Polyhedron('lb',[2.0,-0.1],'ub',[3.0,0.2]) ...
    );

    %% Algorithm
    if nargin > 0
        v_index = 1;
        switch varargin{v_index}
            case 'Theta'
                fs_settings.Theta = varargin{v_index+1};
                v_index = v_index + 2;
            case 'theta'
                fs_settings.theta = varargin{v_index+1};
                if ~fs_settings.Theta.contains(fs_settings.theta)
                    error(['The theta input that was given is not inside the set ''Theta''. Please propose a new one'])
                end
                v_index = v_index + 2;
            case 'm'
                fs_settings.m = varargin{v_index+1};
                v_index = v_index + 2;
            otherwise
                error(['Unexpected input to SimpleSystem1: ' varargin{v_index} ])
        end

    end

    %% Create Outputs
    Theta = fs_settings.Theta;
    U = fs_settings.U;
    X0 = fs_settings.X0;
    X_T = fs_settings.X_T;
end
