function [theta_NMLE, lam] = GASFacCop_G(copula_data, factor_dist, innov_dist, model, max_iter, group_code)
    
    % This code is used to estimate a factor copula with GAS dynamics when the
    % total number of groups G is strictly less than the number of variables N.
    
    % The case G = N is dealt in the file A_main_GASFacCop_N. 
    
    
    Ngroup = max(group_code);
    
    % Nodes and weights for Gauss-Legandre quadrature. 
    % Refer to Oh and Patton (2013) for more details.
    % In order to compute the copula density we need to compute one dimensional
    % integrals through quadrature. We fix the number of nodes and get the
    % corresponding weights.
    Npt_Quad = 5; 
    [Wa, Wb] = GLNodeWt(Npt_Quad);
    GLweight = [Wa, Wb];
    
    %%% Should continue execution from here
    
    
    %%% Factor loadings at t = 1. Here we use the first 65 observations in order to fix initial values.
    corr_ini = corr(copula_data(1:65,:));
    rho_mean = nan(Ngroup,1);
    
    for i = 1:Ngroup
        inx = (group_code == i );
        rho_mean(i,1) = mean(nonzeros(triu(corr_ini(inx,inx),1)));
    end
    
    loading_ini = sqrt(abs(rho_mean'./(1-rho_mean')));
    
    options = optimset('Display','iter','MaxIter',max_iter,'TolFun',1e-6,'TolX',1e-6);
    
    bounds = getBounds(factor_dist, innov_dist, model, Ngroup);
    
    
    %%% Estimating factor copula with GAS recursion
    tic
    
    theta_ini =  [  zeros(Ngroup,1) - 1 ;0.9  ; 0.1 ;  0.1;  0.9; -0.8];
    [theta_NMLE, fobj, exitflag, output1]= fminsearchbnd('LLSum_GASFacCop_G', theta_ini, bounds.lower, bounds.upper, options, copula_data, GLweight, group_code, loading_ini);
    [obj, LL, lam, log_lam, s] = LLSum_GASFacCop_G(theta_NMLE, copula_data, GLweight, group_code, loading_ini);
    
    esti_time = toc  
    
end
