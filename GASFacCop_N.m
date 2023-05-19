function [theta_NMLE, lam] = GASFacCop_N(copula_data, factor_dist, innov_dist, model, max_iter)
    
    % This code is used to estimate a factor copula with GAS dynamics when the
    % total number of groups G is equal to the number of variables N.
    
    %%% This is used only for G = N
    
   % Nodes and weights for Gauss-Legandre quadrature. 
    % Refer to Oh and Patton (2013) for more details.
    % In order to compute the copula density we need to compute one dimensional
    % integrals through quadrature. We fix the number of nodes and get the
    % corresponding weights. 
    Npt_Quad = 50;               
    [Wa, Wb] = GLNodeWt(Npt_Quad);
    GLweight = [Wa, Wb];
    
    %%% Factor loadings at t = 1
    cr_65 = corr(copula_data(1:65,:),'type','spearman');
    cr_full = corr(copula_data,'type','spearman');
    
    lam_ini   = rhobar2betabar(cr_65) ;
    lam_bar = rhobar2betabar(cr_full) ;
    
    %%% Bound for parameters [alpha, beta, nuinv_z, nuinv_eps, psi_z] 
    bounds = getBounds(factor_dist, innov_dist, model, 0);

    options = optimset('Display','iter','MaxIter',max_iter,'TolFun',1e-4,'TolX',1e-4);
    %lower     =  [ 0.0001 ; 0.01;      0.008;  0.008; -0.6 ];
    %upper     =  [ 0.2    ; 0.999999;  0.42;   0.42;   0.6 ]; 
    
    %%% Estimating factor copula with GAS recursion
    tic
    
    theta_ini =  [  0.05  ; 0.95 ;  0.1;  0.1; 0.1];
    [theta_NMLE, ~, ~, ~]= fminsearchbnd('LLSum_GASFacCop_N', theta_ini, bounds.lower, bounds.upper, options, copula_data, GLweight, lam_bar, lam_ini);
    [~, ~, lam, ~, ~] = LLSum_GASFacCop_N(theta_NMLE, copula_data, GLweight, lam_bar, lam_ini);
    
    esti_time = toc

end