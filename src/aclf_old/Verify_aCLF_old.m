function [ Va, dVa_dx_cand, dVa_dth_cand ] = Verify_aCLF_old( sym_x, sym_theta , Va_cand , capa1_sys , Gamma_in , rho0 )
    %VERIFY_ACLF_OLD Verifies aCLFs in the form introduced in Krstic et al. and
    %Taylor et al.
    %   Description
    %   
    %   Usage
    %
    %   Inputs
    %       capa1_sys: This is a control affine parameter affine system of
    %           the form
    %               dx/dt = f(x) + F(x) theta + g(x) u
    %           We require the system to have methods .f(), .F(), and .g()
    %           (reflecting the dynamics) as well as the fields .Theta and 
    %           .U which describe the parameter and input sets
    %           respectively.


    % Constants

    V_Theta = capa1_sys.Theta.V;
    n_VTheta = size(V_Theta,1);

    V_U = capa1_sys.U.V;
    n_VU = size(V_U,1);

    dim_x = length(sym_x);
    dim_theta = length(sym_theta);

    if ~exist('rho0')
        rho0 = 100;
    end

    decay_rate = 0;

    % Create Program
    % ==============
    syms offset0;
    Program1 = sosprogram([sym_x;sym_theta]);

    % Create Constraints
    % ==================

    % Create extra polynomial with adjustable parameter
    offset_poly = -sym_x'*sym_x+offset0;

    % Va should be p.d. in x, but not necessarily in theta

    monom1 = monomials([sym_x],[0,4]);
    [Program1, pd1] = sossosvar(Program1,monom1);

    monom2 = monomials([sym_x;sym_theta],[0:2]);
    [Program1,u1] = sospolyvar(Program1,monom2);

    %Program1 = sosineq( Program1 , Va_cand - pd1 ); % Va_cand - pd1 >= 0

    % Va satisfies gradient condition

    dVa_dx_cand = gradient(Va_cand,sym_x)';
    dVa_dth_cand = gradient(Va_cand,sym_theta)';

    % Create gradient condition
    lambda0_set = {};
    set_of_Lambdas = {};

    % Compute the gradient condition for each theta.
    [Program1, lambda0] = sossosvar(Program1,monom1); %reate lambda 1
    lambda0_set{1} = lambda0;

    % Start compiling expression
    grad_like_term = (1+lambda0)*(Va_cand - rho0)*(sym_x'*sym_x);

    for theta_index = [1:n_VTheta]

        % Get the theta at this vertex
        v_Theta = V_Theta(theta_index,:)';

        Lambda = {};
        for v_U_Index = [1:n_VU]
            [Program1, Lambda{v_U_Index}] = sossosvar(Program1,monom1);
        end
        set_of_Lambdas{theta_index} = Lambda;
        
        for v_U_index = [1:n_VU]
            v_U = V_U(v_U_index,:)';

            % Add constraints for every input
            dVa_dt_cand = dVa_dx_cand * ( ...
                capa1_sys.f(sym_x) + ...
                capa1_sys.F(sym_x) * (v_Theta + Gamma_in * dVa_dth_cand' ) + ...
                capa1_sys.g(sym_x) * v_U ...
                );

            for u_index = 1:n_VU
                grad_like_term = grad_like_term + ...
                    Lambda{v_U_index}*(dVa_dt_cand + decay_rate * Va_cand);
            end

        end

%         % Enforce gradient is always decreasing.
%         Program1 = sosineq( Program1, -grad_like_term + offset_poly ) % -grad_like_term >= 0  -> grad_like_term <= 0

    end

    % Enforce gradient is always decreasing.
    Program1 = sosineq( Program1, grad_like_term ) % -grad_like_term >= 0  -> grad_like_term <= 0

    % Solve
    % =====

    Program1 = sossolve( Program1 );
    assert( (Program1.solinfo.info.pinf ~= 1) || (Program1.solinfo.info.dinf ~= 1) ); % Make sure that neither problem is infeasible?

    Va = sosgetsol(Program1,Va_cand); %Getting solution for V
    dVa_dx_cand = sosgetsol(Program1,dVa_dx_cand); %Getting solution for dha_dx
    dVa_dth_cand = sosgetsol(Program1,dVa_dth_cand); %Getting solution for dha_dx

end

