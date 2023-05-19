function out1 = X_cdf(x, theta, GLweight)

% The marginal CDF of X_i associated with skew t - t factor model
%
% INPUTS:   x, a scalar, the value of x that we will evaluate G at
%           theta, =[lam, nuinv_z, nuinv_esp, psi_z], the parameters of Fz and Feps
%           GLweight, a matrix, nodes and weights for Gauss-Legendre quadrature
%
%  OUTPUTS: out1, a scalar, the value of the marginal CDF at x
 

out1 = GLquad('X_cdf_helper', 1e-5, 1-1e-5, GLweight, x, theta);

end