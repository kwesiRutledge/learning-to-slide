classdef SOSTest < matlab.unittest.TestCase

    methods (TestClassSetup)
        % Shared setup for the entire test class
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        % Test methods

        function yalmipTutorialTest1(testCase)
            % Set up the basic SOS optimization from the YALMIP Tutorial
            x = sdpvar(1,1);y = sdpvar(1,1);
            p = (1+x)^4 + (1-y)^2;
            F = sos(p);
            solvesos(F);

            testCase.verifyTrue(true);
        end
    end

end