classdef SOSToolsTest < matlab.unittest.TestCase
    
    methods(TestClassSetup)
        % Shared setup for the entire test class
    end
    
    methods(TestMethodSetup)
        % Setup for each test
    end
    
    methods(Test)
        % Test methods
        
        function basicSOSVariableTest(testCase)
            % Create a dummy program from the tutorial for SOSTools
            syms x y a b c d;
        
            Program1 = sosprogram([x;y]);
            [Program1,p1] = sossosvar(Program1,[x;y]);
            
            %p = 2*x^2 + 3 * x * y + 4 * y^4;
            p = 2*x^2 + 4 * y^4;

            %Add Constraints
            Program1 = sosineq( Program1 , p );

            Program1 = sossolve(Program1)
            Program1.solinfo.info.pinf

            testCase.verifyTrue( ...
                (Program1.solinfo.info.pinf == 0) && (Program1.solinfo.info.dinf == 0) ...
            );
        end

        function yalmipSOSVariableTest(testCase)
            % Create a dummy program from the tutorial for YALMIP
            syms x y a b c d;
        
            Program1 = sosprogram([x;y]);
            [Program1,p1] = sossosvar(Program1,[x;y]);
            
            %p = 2*x^2 + 3 * x * y + 4 * y^4;
            p = (1+x)*x^4 + (1-y)^2;

            %Add Constraints
            Program1 = sosineq( Program1 , p );
            Program1 = sossolve( Program1 );

            testCase.verifyTrue( ...
                (Program1.solinfo.info.pinf ~= 0) || (Program1.solinfo.info.dinf ~= 0) ...
            );
        end

        function yalmipSOSVariableTest2(testCase)
            % Constants
            Theta1 = Polyhedron('lb',0.5,'ub',0.85);
            system1 = SimpleSystem1(Theta1)

            % Create a dummy program from the tutorial for YALMIP
            syms x a b c d;
        
            Program1 = sosprogram([x]);
            monom1 = monomials([x],[1,2,3,4]);
            [Program1,V] = sospolyvar(Program1,monom1);
            
            l1 = 0.01*x^2;
            l2 = 0.01*x^2;
            
            s1 = x^2;
            [Program1,p2] = sospolyvar(Program1,monom1);
            
            dV_dx = diff(V,x);
            
            % Create constraints
            Program1 = sosineq( Program1 , V - l1 ); %F2
            tempExpr = ...
                - ( s1 * dV_dx * ( system1.f(x) + system1.F(x)* system1.theta ) + ...
                p2 * (dV_dx*( system1.g(x) + system1.G(x) * system1.theta ) ) + l2)
            Program1 = sosineq( Program1 , ...
                tempExpr ...
            ); %F3

            % Optimize
            Program1 = sossolve( Program1 );

            testCase.verifyTrue( ...
                (Program1.solinfo.info.pinf ~= 0) || (Program1.solinfo.info.dinf ~= 0) ...
            );
        end
        
    end
    
end