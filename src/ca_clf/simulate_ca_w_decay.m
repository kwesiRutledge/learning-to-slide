function [x,u_history, V_history] = simulate_ca_w_decay(ca_sys, V_symb, T_sim, X0 , dt, N_sims)
%SIMULATE_CA_W_DECAY Summary of this function goes here
%   Detailed explanation goes here

    % Constants
    decay_rate = 1;

    T_sim = 10; %Number of steps to simulate system forward.
    dt = 0.1;
    
    if ~exist('N_sims')
        N_sims = 1;
    end

    x0 = sampleFromPolytope(X0,N_sims);
    n_x = X0.Dim;

    % Compute gradient of V_symb

    dV_dx_symb = gradient(V_symb)';

    % Simulate

    x = {}; u_history = {}; V_history = {};
    for sim_index = 1:N_sims
        x_i = [x0(:,sim_index)];
        u_i = [];
        V_i = double(subs(V_symb,struct('x1', x_i(1), 'x2', x_i(2))));
    
        for k = [0:T_sim-1]
            xi_k = x_i(:,k+1);
            
            Vi_k = double(subs(V_symb,struct('x1', xi_k(1), 'x2', xi_k(2))));
            dVi_k_dx = double(subs(dV_dx_symb,struct('x1', xi_k(1), 'x2', xi_k(2))));
    
            [ ui_k, optim_out ] = CLF_control_ca_systems(xi_k , Vi_k, dVi_k_dx, ca_sys, decay_rate);
            xi_kp1 = xi_k + ca_sys.dynamics(xi_k,ui_k)*dt;
    
            Vi_kp1 = double(subs(V_symb,struct('x1', xi_kp1(1), 'x2', xi_kp1(2))));
    
            % Update the tracking
            x_i = [x_i, xi_kp1];
            u_i = [u_i, ui_k];
            V_i = [V_i, Vi_kp1];
    
        end
    
        x{sim_index} = x_i;
        u_history{sim_index} = u_i;
        V_history{sim_index} = V_i;
    
    end

end

