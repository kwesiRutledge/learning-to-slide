function [x, th0, V_history, th_history, th_star, u_history] = simulate_capa2(capa2_sys, Va_symb, T_sim, X0 , dt, Gamma, N_sims, decay_rate, state_names, theta_names)
    %SIMULATE_CAPA2 Simulates the control parameter affine system
    %(version 2) with the aCLF Va_symb
    %Usage
    %
    %Inputs
    %   capa2_sys: A control affine parameter affine system of the form
    %           dxdt = f(x) + F(x) * theta + g(x) u
    %       It is assumed that the system has the following methods:
    %           f(), F(), g(), dynamics()
    %       and the following fields
    %           Theta, theta
    
    %% Constants
    x_dim = X0.Dim;
    theta_dim = capa2_sys.Theta.Dim;

    %% Input Processing
    if ~exist('N_sims')
        N_sims = 1;
    end

    if ~exist('decay_rate')
        decay_rate = 0;
    end

    if ~exist('state_names')
        state_names = {};
        for x_index = 1:x_dim
            state_names{x_index} = ['x' num2str(x_index)];
        end
    end

    if ~exist('theta_names')
        theta_names = {};
        for theta_index = 1:theta_dim
            theta_names{theta_index} = ['theta' num2str(x_index)];
        end
    end
    
    %% Compute Derivatives
    for x_index = 1:x_dim
        x_symb(x_index) = sym(state_names(x_index));
    end
    dVa_dx_symb = gradient(Va_symb,x_symb)';
    
    for theta_index = 1:theta_dim
        theta_symb(theta_index) = sym(theta_names(theta_index));
    end
    dVa_dth_symb = gradient(Va_symb,theta_symb)';

    %% Simulate
    
    % Sample
    x0 = sampleFromPolytope(X0,N_sims);
    th0 = sampleFromPolytope(capa2_sys.Theta, N_sims);
    th_star = sampleFromPolytope(capa2_sys.Theta, N_sims);
    
    x = {}; u_history = {}; V_history = {}; th_history = {};
    for sim_index = 1:N_sims
        %Set Initial Conditions
        x_i = [x0(:,sim_index)];
        th_i = [th0(:,sim_index)];
        capa2_sys.theta = th_star(:,sim_index); %th_i;

        u_i = [];
        V_i = double( ...
            subs( ...
                Va_symb, ...
                state_parameter_pair_to_struct(x_i(:,1),th_i(:,1),state_names,theta_names) ...
            ));
    
        for k = [0:T_sim-1]
            xi_k = x_i(:,k+1);
            thi_k = th_i(:,k+1);
    
            Vi_k = double( ...
                subs( ...
                    Va_symb, ...
                    state_parameter_pair_to_struct(xi_k,thi_k,state_names,theta_names) ...
                ));
            dVi_k_dx = double( ...
                subs( ...
                    dVa_dx_symb, ...
                    state_parameter_pair_to_struct(xi_k,thi_k,state_names,theta_names) ...
                ));
    
            [ ui_k, optim_out ] = aCLF_control_capa2( ...
                xi_k , ...
                Vi_k, dVi_k_dx, dVa_dth_symb, ...
                capa2_sys, ...
                decay_rate, Gamma, ...
                state_names, theta_names);

            if optim_out.problem ~= 0
                warning('Warning! CLF control optimization was not feasible!')
                optim_out;
            end

            xi_kp1 = xi_k + capa2_sys.dynamics(xi_k,ui_k)*dt;
            thi_kp1 = thi_k + aCLF_estimation_capa2(xi_k,ui_k,thi_k,dVi_k_dx,capa2_sys)*dt; 
    
            Vi_kp1 = double( ...
                subs( ...
                    Va_symb, ...
                    state_parameter_pair_to_struct(xi_kp1,thi_kp1,state_names,theta_names) ...
                ));
    
            % Update the tracking
            x_i = [x_i, xi_kp1];
            th_i = [th_i, thi_kp1];
            u_i = [u_i, ui_k];
            V_i = [V_i, Vi_kp1];
    
        end
    
        x{sim_index} = x_i;
        u_history{sim_index} = u_i;
        V_history{sim_index} = V_i;
        th_history{sim_index} = th_i;
    
    end

end