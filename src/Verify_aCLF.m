function [ Va, dVa_dx, dVa_dth ] = Verify_aCLF( sym_x, sym_theta , Va_cand , cada_sys , Gamma_in )
    %VERIFY_ACLF Verifies the input symbolic aCLF.
    %Description:
    %   Tests to make sure that the provided aCLF (symbolic function)
    %   satisfies the condition of an aCLF for control-affine,
    %   parameter-affine systems:
    %       dxdt = f(x) + F(x) \theta + (g(x) + \sum_i \theta_i G_i(x)) u
    %
    %Usage
    %
    %Inputs
    %   V_cand: A symbolic expression which should be a candidate aCLF.
    %   cada_sys: A Control-Affine, Parameter-Affine System defined with
    %               the assumptions I've baked into some examples like
    %               scalar_system.
    %   sym_x: A vector of symbolic containing one variable for each state.

    % Constants

    V_Theta = cada_sys.Theta.V;
    n_VTheta = size(V_Theta,1);

    dim_x = length(sym_x);
    dim_theta = length(sym_theta);

    % Create Program
    % ==============
    syms offset0;
    Program1 = sosprogram([sym_x;sym_theta],[offset0]);

    % Create Constraints
    % ==================

    % Create extra polynomial with adjustable parameter
    offset_poly = -sym_x'*sym_x+offset0;

    % Va should be p.d. in x, but not necessarily in theta

    monom1 = monomials([sym_x],[0,4]);
    [Program1, pd1] = sossosvar(Program1,monom1);

    monom2 = monomials([sym_x;sym_theta],[0:2]);
    [Program1,u1] = sospolyvar(Program1,monom2);

    Program1 = sosineq( Program1 , Va_cand - pd1 ); % Va_cand - pd1 >= 0

    % Va satisfies gradient condition

    for x_index = [1:dim_x]
        dVa_dx_cand(x_index,1) = diff(Va_cand,sym_x(x_index));
    end
    for theta_index = [1:dim_theta]
        dVa_dth_cand(theta_index,1) = diff(Va_cand,sym_theta);
    end

    u1 = -1*(sym_x-(-2.5*sym_theta - (13/4)) );

    for v_index = [1:n_VTheta]
        % Get the theta at this vertex
        v_Theta = V_Theta(v_index,:)';

        % Create Gsum term
        G = cada_sys.G(sym_x);
        sum_Gi = v_Theta(1)*G{1};
        for theta_dim = [ 2 : dim_theta ]
            sum_Gi = sum_Gi + v_Theta(theta_dim)*G{theta_dim};
        end

        grad_like_term = dVa_dx_cand' * ( ...
            cada_sys.f(sym_x) + cada_sys.F(sym_x) *( v_Theta + Gamma_in * dVa_dth_cand) + ...
            (cada_sys.g(sym_x) + sum_Gi)*u1 ...
            );

        % Enforce gradient is always decreasing.
        Program1 = sosineq( Program1, -grad_like_term + offset_poly ) % -grad_like_term >= 0  -> grad_like_term <= 0

    end

    % Solve
    % =====

    Program1 = sossolve( Program1 );
    assert( (Program1.solinfo.info.pinf ~= 1) || (Program1.solinfo.info.dinf ~= 1) ); % Make sure that neither problem is infeasible?

    Va = sosgetsol(Program1,Va_cand); %Getting solution for V
    dVa_dx = sosgetsol(Program1,dVa_dx_cand); %Getting solution for dha_dx
    dVa_dth = sosgetsol(Program1,dVa_dth_cand); %Getting solution for dha_dx

end