classdef AdaptiveCruiseControl
    %ADAPTIVECRUISECONTROL Simple ACC System (2-Dimensional)
    %   Detailed explanation goes here
    
    properties
        v0
        theta
        m
    end
    
    methods
        function obj = AdaptiveCruiseControl()
            %ADAPTIVECRUISECONTROL Construct an instance of this class
            %   I don't think we will allow people to set these values.
            obj.v0 = 20; %20 m/s
            obj.theta = [   51;     % N
                            1.2567; % Ns/m
                            0.4342];% Ns^2/m^2
            obj.m = 1370; % kg
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
            v0 = obj.v0;

            % Algorithm
            v = x(1);
            f_out = [ 0 ; v0 - v ];
        end

        function F_out = F(obj,x)
            %F State-dependent Linear Term on dynamics equation
            %Description:
            %   Computes the F() term in the dynamics.

            % Constants
            v0 = obj.v0;
            m = obj.m;

            % Algorithm
            v = x(1);
            F_out = (1/m)*[ 1, v, v.^2 ; 0 , 0 , 0 ];
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

