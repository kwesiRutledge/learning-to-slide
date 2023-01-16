function [ u_opt , optim_out ] = aCLF_control_capa2( x, Va, dVa_dx, dVa_dth_symb, capa2_sys, decay_rate, Gamma_in, state_names, theta_names )
    %CLF_CONTROL_CAPA2 Computes aCLF control of the system
    %Description
    %
    %Usage
    %   [ u_opt , optim_out ] = aCLF_control_capa2( ...
    %                               x, Va, dVa_dx, dVa_dth_symb, capa2_sys, ...
    %                               'decay_rate' , decay_rate , ...
    %                               'Gamma', Gamma_in , ...
    %                               'state_names', state_names, ...
    %                               'theta_names', theta_names)
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
    %               Theta, U, dim_theta, dim_x, dim_u

    % Input Processing
    %[x, Va, dVa_dx, dVa_dth_symb, capa2_sys, control_settings] = input_processing_aCLF_control_capa2(varargin{:})
    
%     decay_rate = control_settings.decay_rate;
% 
%     Gamma = control_settings.Gamma;
% 
%     state_names = control_settings.state_names;
%     theta_names = control_settings.theta_names;

    if ~exist("decay_rate")
        decay_rate = 0;
    end

    if ~exist("Gamma_in")
        Gamma_in = eye(capa2_sys.dim_theta);
    end

    if ~exist("theta_names")
       theta_names = cell(1,capa2_sys.dim_theta);
        for theta_index = 1:capa2_sys.dim_theta
            theta_names{theta_index} = ['theta' num2str(theta_index)];
        end
    end

    if ~exist("state_names")
        state_names = cell(1,capa2_sys.dim_x);
        for x_index = 1:capa2_sys.dim_x
            state_names{x_index} = ['x' num2str(x_index)];
        end
    end


    % Constants
    settings_in = sdpsettings('verbose',0);
    dim_x = capa2_sys.dim_x;
    dim_u = capa2_sys.dim_u;

    Theta = capa2_sys.Theta;
    dim_theta = Theta.Dim;
    V_Theta = Theta.V;
    n_VTheta = size(V_Theta,1);

    U = capa2_sys.U;

    % Algorithm
    % =========

    u = sdpvar(dim_u,1);
    constraints = [];

    Lf_Va = dVa_dx * capa2_sys.f(x);
    Lg_Va = dVa_dx * capa2_sys.g(x);
    LF_Va = dVa_dx * capa2_sys.F(x);
    
    LG_Va = {};
    G = capa2_sys.G(x);
    for theta_index = 1:dim_theta
        LG_Va{theta_index} = dVa_dx * G{theta_index};
    end

    % Dynamics Constraints
    for theta_index = 1:n_VTheta

        v_Theta = V_Theta(theta_index,:);
        
        dVa_dth_xv = double( ...
            subs( ...
                dVa_dth_symb, ...
                state_parameter_pair_to_struct(x,v_Theta,state_names,theta_names) ...
            ));

        sum_Gi = v_Theta(1)*LG_Va{1};
        for theta_index = 2:dim_theta
            sum_Gi = sum_Gi + v_Theta(theta_index) * LG_Va{theta_index};
        end

        lyapunov_function_decrease_constraint = ...
            [ Lf_Va + LF_Va * (v_Theta + Gamma_in*dVa_dth_xv') + (Lg_Va + sum_Gi)*u <= - decay_rate * Va];

        constraints = constraints + lyapunov_function_decrease_constraint;

    end

    % Input Bound Constraints
    input_constraint = [capa2_sys.U.A * u <= capa2_sys.U.b];
    constraints = constraints + input_constraint;

    % Optimize
    optim_out = optimize(constraints, u'*u ,settings_in);

    u_opt = value(u); % Save the feasible u

end

function [x, Va, dVa_dx, dVa_dth_symb, capa2_sys, control_settings] = input_processing_aCLF_control_capa2(varargin)
    %Description:
    %   Process the inputs given to ExternalBehaviorSet() constructor.

    %% Required arguments
    if nargin < 5
        error(['Require at least 5 arguments; received ' num2str(nargin) ])
    end

    x = varargin{1};
    Va = varargin{2};
    dVa_dx = varargin{3};
    dVa_dth_symb = varargin{4};
    capa2_sys = varargin{5};

    %% Set Defaults
    state_names = cell(1,capa2_sys.dim_x);
    for x_index = 1:capa2_sys.dim_x
        state_names{x_index} = ['x' num2str(x_index)];
    end
    theta_names = cell(1,capa2_sys.dim_theta);
    for theta_index = 1:capa2_sys.dim_theta
        theta_names{theta_index} = ['theta' num2str(theta_index)];
    end

    control_settings = struct( ...
        'decay_rate', 0, ...
        'Gamma', eye(capa2_sys.dim_theta), ...
        'state_names', state_names, ...
        'theta_names', theta_names ...
    );

    %% Algorithm
    if nargin > 0
        v_index = 6;
        switch varargin{v_index}
            case 'theta_names'
                control_settings.theta_names = varargin{v_index+1};
                if length(control_settings.theta_names) ~= capa2_sys.dim_theta
                    error(['Expected number of state names (' num2str(length(control_settings.theta_names)) ') and dimension of unknown parameters (' num2str(capa2_sys.dim_theta) ') to be the same. They aren''t!' ])
                end
                v_index = v_index + 2;
            case 'state_names'
                control_settings.state_names = varargin{v_index+1};
                if length(control_settings.state_names) ~= capa2_sys.dim_x
                    error(['Expected number of state names (' num2str(length(control_settings.state_names)) ') and dimension of system (' num2str(capa2_sys.dim_x) ') to be the same. They aren''t!' ])
                end
                v_index = v_index + 2;
            case 'decay_rate'
                control_settings.decay_rate = varargin{v_index+1};
                v_index = v_index + 2;
            otherwise
                error(['Unexpected input to SimpleSystem1: ' varargin{v_index} ])
        end

    end

    %% Create Outputs
    
end