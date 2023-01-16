function [ tau ] = aCLF_estimation_capa2(x,u,theta_hat, dVa_dx,capa2_sys)
    %ACLF_ESTIMATION_CAPA2 This is the tau function governing the
    %dynamics of the estimator of the unknown parameters theta in capa1
    %systems.
    %   Description
    %
    %   Usage
    %
    %   Inputs
    %   dVa_dx : Current value of the adaptive CLF's gradient with respect
    %           to the state x
    %           Should be a 1 x n vector where n is the dimension of the
    %           state
    %   capa2_sys: Control-Affine Parameter Affine System Type 1
    %                   dxdt = f(x) + F(x) \theta + g(x) u
    %              These systems should have the following methods defined
    %                   f(), F(), g()
    %               and the following member values:
    %                   Theta, U

    %% Input Processing

    dim_theta = capa2_sys.dim_theta;

    G = capa2_sys.G(x);
    for theta_index = 1:dim_theta
        LG_Va{theta_index} = dVa_dx * G{theta_index};
    end

    sum_LG_Vi = LG_Va{1} * theta_hat(1);
    for theta_index = 2:dim_theta
        sum_LG_Vi = sum_LG_Vi + LG_Va{theta_index} * theta_hat(theta_index);
    end

    tau = (dVa_dx * capa2_sys.F(x) + sum_LG_Vi*u)';

end

