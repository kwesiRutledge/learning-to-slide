function [ tau ] = aCLF_estimation_capa1_systems(x,dVa_dx,capa1_sys)
    %ACLF_ESTIMATION_CAPA1_SYSTEMS This is the tau function governing the
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
    %   capa1_sys: Control-Affine Parameter Affine System Type 1
    %                   dxdt = f(x) + F(x) \theta + g(x) u
    %              These systems should have the following methods defined
    %                   f(), F(), g()
    %               and the following member values:
    %                   Theta, U

    %% Input Processing

    tau = (dVa_dx * capa1_sys.F(x))';

end

