function [ u_opt, optim_out ] = CLF_control_ca_systems( x, V_x, dV_dx, ca_sys, decay_rate )
    %CLF_CONTROL_CA_SYSTEMS Summary of this function goes here
    %Description
    %
    %Inputs
    %   V : Current value of the adaptive CLF
    %   dV_dx : Current value of the adaptive CLF's gradient with respect
    %           to the state x
    %           Should be a 1 x n vector where n is the dimension of the
    %           state
    %   


    % Constants
    settings_in = sdpsettings('verbose',0);
    dim_x = ca_sys.dim_x;
    dim_u = ca_sys.dim_u;

    if ~exist('decay_rate')
        decay_rate = 0;
    end

    % Algorithm
    % =========

    u = sdpvar(dim_u,1);
    constraints = [];

    Lf_V = dV_dx * ca_sys.f(x);
    Lg_V = dV_dx * ca_sys.g(x);
    
    % Dynamics Constraint
    lyapunov_function_decrease_constraint = [Lf_V + Lg_V*u <= - decay_rate * V_x];
    
    constraints = constraints + lyapunov_function_decrease_constraint;

    % Input Bound Constraints
    input_constraint = [ca_sys.U.A * u <= ca_sys.U.b];
    constraints = constraints + input_constraint;

    % Optimize
    optim_out = optimize(constraints, u'*u ,settings_in);

    u_opt = value(u); % Save the feasible u

end

