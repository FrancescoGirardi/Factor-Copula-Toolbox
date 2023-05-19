

%%% INPUT: Distributions for the factor and innovations.
%%% OUTPUT: lower and upper bounds for parameters.


function [out] = getBounds(factor_dist, innov_dist, model, Ngroup)

    %%Input check:
    
    factor_check = strcmp(factor_dist,'Normal') | strcmp(factor_dist,'t') | strcmp(factor_dist,'skew_t');
    if factor_check == false
        error('Invalid distribution for factor.')
    end
    
    innov_check = strcmp(innov_dist,'Normal') | strcmp(innov_dist,'t');
    if innov_check == false
        error('Invalid distribution for innovation.')
    end

    model_check = strcmp(model,'Dynamic') | strcmp(model,'Static');
    if model_check == false
        error('Invalid model type')
    end

   
    % Bound for parameters [omega's, alpha, beta, nuinv_z, nuinv_eps, psi_z] 
    % nuinv_z: Inverse of degrees of freedom for factor.
    % nuinv_eps: Inverse of degrees of freedom for innovation.
    % psi_z: Non-centrality parameter for factor.

    %%% Standard parameters for full generality.
    lower     =  [ -0.3*ones(Ngroup,1); 0.0001 ; 0.01;      0.008;  0.008; -0.6 ];
    upper     =  [  0.3*ones(Ngroup,1); 0.2    ; 0.999999;  0.42;   0.42;   0.6 ]; 
    
    if strcmp(factor_dist,'Normal')
        %Set degrees of freedom to infinity and non centrality to 0
        lower(end - 2) = 0.001;
        upper(end - 2) = 0.01;
        lower(end) = 0;
        upper(end) = 0;
    end

    if strcmp(factor_dist,'t')
        %Set non centrality parameter to 0
        lower(end) = 0;
        upper(end) = 0;
    end

    if strcmp(innov_dist,'Normal')
        %Set degrees of freedom to infinity
        lower(end - 1) = 0.05;
        upper(end - 1) = 0.10;
    end

    if strcmp(model,'Static')
        %Set alpha e beta to 0
        lower(end - 4) = 0;
        upper(end - 4) = 0;
        lower(end - 3) = 0;
        upper(end - 3) = 0;
    end

    out.lower = lower;
    out.upper = upper;

end