%% System Definition
classdef DaiP_ToySystem1
    properties
        U
        dim_x
        dim_u
    end

    methods
        function [obj] = DaiP_ToySystem1()
            %Description
            %   Initializes the two dimensional toy system.
    
            obj.U = Polyhedron('lb',-0.4,'ub',0.4);
            obj.dim_x = 2;
            obj.dim_u = 1;
    
        end
    
        function [f_out] = f(obj,x)
            %Description
            %
            %Usage
            %   f_out = ts.f(x)
    
            f_out = [0;-x(1) + (1/6)*x(1)^3];
        end
    
        function [g_out] = g(obj,x)
            %Description
            %
            %Usage
            %   g_out = ts.g(x)
    
            g_out = [1;-1];
        end
    
        function [x_dot] = dynamics(obj,x,u)
            %Description
            %   Computes the derivative of the state w.r.t. time given the
            %   current state and input u.
            
            x_dot = obj.f(x) + obj.g(x)*u;
    
        end
    end

end