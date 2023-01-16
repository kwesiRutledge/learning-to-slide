classdef AdaptiveCruiseControlwNominal
    %ADAPTIVECRUISECONTROLWNOMINAL Simple ACC System (2-Dimensional) with
    %target velocity and headway value baked into the state.
    %   Description
    %       Dynamics are:
    %           \dot{x} = [   0   ] + (1/m) [ 1, (v-v*) , (v-v*)^2 ] [ f1 ] + [ 1/m ] u 
    %                     [ v0 - v]         [ 0,   0    ,    0     ] [ f2 ] + [  0  ] 
    %                                                                [ f3 ]
    
    properties
        % System Known Parameters
        v0;
        m;
        U;
        dim_x;
        dim_u;
        D_star;
        v_star;
        % System's Unknown Parameters
        theta;
        Theta;
    end
    
    methods
        function obj = AdaptiveCruiseControlwNominal()
            %ADAPTIVECRUISECONTROL Construct an instance of this class
            %   I don't think we will allow people to set these values.
            obj.v0 = 20; %20 m/s
            obj.theta = [   51;     % N
                            1.2567; % Ns/m
                            0.4342];% Ns^2/m^2
            obj.m = 1370; % kg
            
            dtheta = [1.0,0.1,0.02];
            obj.Theta = Polyhedron('lb',-dtheta,'ub',dtheta) + obj.theta;

            obj.U = Polyhedron('lb',-2*9.8*obj.m,'ub',2*9.8*obj.m);

            % Dimension
            obj.dim_x = 2;
            obj.dim_u = 1;

            % Nominal states
            obj.v_star = obj.v0;
            obj.D_star = 20;
        end
        
        function dx_dt = dynamics(obj,x,u)
            %dynamics
            %Description:
            %   The dynamics function constructs the derivative of the
            %   state which is:
            %       dx_dt = f(x) + F(x) * theta + g(x) u

            dx_dt = obj.f(x) + obj.F(x) * obj.theta + obj.g(x)*u;
        end

        function f_out = f(obj,x)
            %f State-dependent offset on dynamics equation
            %Description:
            %   Computes the f() term in the dynamics.

            % Constants

            % Algorithm
            v_bar = x(1);
            f_out = [ 0 ; obj.v0 - v_bar - obj.v_star ];
        end

        function F_out = F(obj,x)
            %F State-dependent Linear Term on dynamics equation
            %Description:
            %   Computes the F() term in the dynamics.

            % Constants
            m = obj.m;

            % Algorithm
            v_bar = x(1);
            F_out = -(1/m)*[ 1, (v_bar+obj.v_star), (v_bar+obj.v_star).^2 ; 0 , 0 , 0 ];
        end

        function g_out = g(obj,x)
            %g State-dependent control gain in dynamics equation
            %Description:
            %   Computes the g() term in the dynamics.

            % Constants
            m = obj.m;

            % Algorithm
            g_out = [ 1/m ; 0 ];
        end
    end
end

