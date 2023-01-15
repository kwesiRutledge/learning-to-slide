function [ u_opt, optim_out ] = CLF_control_capa1_systems( x, Va, dVa_dx, dVa_dth_symb, capa1_sys, decay_rate , Gamma_in )
    %CLF_CONTROL_CAPA1_SYSTEMS Summary of this function goes here
    %Description
    %
    %Inputs
    %   V : Current value of the adaptive CLF
    %   dV_dx : Current value of the adaptive CLF's gradient with respect
    %           to the state x
    %           Should be a 1 x n vector where n is the dimension of the
    %           state
    %   capa1: Control-Affine Parameter Affine System Type 1
    %               dxdt = f(x) + F(x) \theta + g(x) u
    %           These systems should have the following methods defined
    %               f(), F(), g()
    %           and the following member values:
    %               Theta, U


    % Constants
    settings_in = sdpsettings('verbose',0);
    dim_x = capa1_sys.dim_x;
    dim_u = capa1_sys.dim_u;

    if ~exist('decay_rate')
        decay_rate = 0;
    end

    Theta = capa1_sys.Theta;
    V_Theta = Theta.V;
    n_VTheta = size(V_Theta,1);

    U = capa1_sys.U;


    % Algorithm
    % =========

    u = sdpvar(dim_u,1);
    constraints = [];

    Lf_V = dVa_dx * capa1_sys.f(x);
    Lg_V = dVa_dx * capa1_sys.g(x);
    LF_V = dVa_dx * capa1_sys.F(x);
    
    % Dynamics Constraints
    for theta_index = 1:n_VTheta

        v_Theta = V_Theta(theta_index,:);
        
        dVa_dth_xv = double(subs(dVa_dth_symb, state_parameter_pair_to_struct(x,v_Theta)));

        lyapunov_function_decrease_constraint = ...
            [ Lf_V + LF_V * (v_Theta + Gamma_in*dVa_dth_xv') + Lg_V*u <= - decay_rate * Va];

        constraints = constraints + lyapunov_function_decrease_constraint;

    end

    % Input Bound Constraints
    input_constraint = [capa1_sys.U.A * u <= capa1_sys.U.b];
    constraints = constraints + input_constraint;

    % Optimize
    optim_out = optimize(constraints, u'*u ,settings_in);

    u_opt = value(u); % Save the feasible u

end

