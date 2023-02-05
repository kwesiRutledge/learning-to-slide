function [V, dV_dx, lambda_set] = Verify_ca_CLF(sym_x, V_cand , ca_sys)
    %VERIFY_CA_CLF Verifies whether or not the candidate Lyapunov function
    %   is a certificate of the origin being a stable equilibrium point
    %   Description
    %   
    %   Inputs
    %       sym_x: A vector of states
    %
    %   Outputs
    %   

    % Constants
    % =========

    dim_x = length(sym_x);
    V_degree = 4;

    V_U = ca_sys.U.V;
    n_VU = size(V_U,1);

    rho0 = 0.1;

    % Create Program
    % ==============
    
    syms offset0 real;
    Program1 = sosprogram([sym_x]);
    
    % Create Constraints
    % ==================
    
    % Create extra polynomial with adjustable parameter
    offset_poly = -sym_x'*sym_x+offset0;
    
    % V should be p.d. in x
    monom1 = monomials([sym_x],[0,V_degree]);
    [Program1, pd1] = sossosvar(Program1,monom1); % Because variable is 
                                                    % defined as SOS, it is
                                                    % positive.
    
    % Va satisfies gradient condition
    
%     for x_index = [1:dim_x]
%         dV_dx_cand(x_index,1) = diff(V_cand,sym_x(x_index));
%     end
    dV_dx_cand = gradient(V_cand, sym_x);
    
    [Program1, lambda0] = sossosvar(Program1,monom1); %reate lambda 1
    lambda_set_cand = {};
    for v_U_Index = [1:n_VU]
        [Program1, lambda_set_cand{v_U_Index}] = sossosvar(Program1,monom1);
    end
    
    % Create gradient condition
    grad_like_expression = (1+lambda0)*(V_cand - rho0)*(sym_x'*sym_x);
    for v_index = [1:n_VU]
        % Get the theta at this vertex
        v_U = V_U(v_index,:)';
    
        grad_like_expression = grad_like_expression - ...
            lambda_set_cand{v_index} * ( dV_dx_cand'*(ca_sys.f(sym_x) + ca_sys.g(sym_x)*v_U) );
    
    end
    % Enforce gradient is always decreasing.
    Program1 = sosineq( Program1, grad_like_expression ) % grad_like_term >= 0 
    
    % Solve
    % =====
    
    Program1 = sossolve( Program1 );
    assert( (Program1.solinfo.info.pinf ~= 1) || (Program1.solinfo.info.dinf ~= 1) ); % Make sure that neither problem is infeasible?
    
    V = sosgetsol(Program1,V_cand); %Getting solution for V
    dV_dx = sosgetsol(Program1,dV_dx_cand); %Getting solution for dha_dx
    lambda_set = {};
    for lambda_index = 1:length(lambda_set_cand)
        lambda_set{lambda_index} = sosgetsol(Program1, lambda_set_cand{lambda_index});
    end

end

